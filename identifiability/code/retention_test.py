"""Test of the two mechanisms of false support on real traditions, for every witness pair (i, j):
  retention ratio   = rate of units where i and j each have a reading of their own and all other legible
                      witnesses share one reading (these propose the unrooted split {i,j}|rest with no shared
                      error), divided by p_i p_j;
  coincidence ratio = rate of units where exactly i and j share a reading against at least two others, divided by
                      p_i p_j (this is kappa-hat of graph_stats.py, here over all pairs);
rates on the units legible in both; p_i = share of units legible in i where i has a reading of its own. Under the
attraction model with independent corruption, the retention ratio is near the probability that the other witnesses
agree (near 1 when corruption is light), whatever alpha, and the coincidence ratio near alpha^2 for non-sister pairs
(true sister pairs share departures by descent and show large coincidence ratios).
Usage: retention_test.py OUT.json"""
import sys, os, json, itertools, statistics as st

HOME = os.path.expanduser('~'); P = os.path.join(HOME, 'MC/projects/priscus_excerpta/data')
SRC = {'Heinrichi': os.path.join(HOME, 'stemma_xfer/heinrichi/heinrichi.num'),
       'Donatus': os.path.join(HOME, 'stemma_xfer/donatus_final.num'),
       'Lawrence (BHL 4759)': os.path.join(P, 'lawrence_bhl4759.num'),
       'Jordanes, Getica': os.path.join(P, 'historic/mommsen_jordanes/jordanes.num'),
       'Iwein (Lachmann 1827)': os.path.join(P, 'historic/iwein_1827/iwein.num')}


def load_num(path):
    L = [l.split() for l in open(path).read().split('\n') if l.strip()]
    return L[0], [r for r in L[1:] if len(r) == len(L[0])]


def load_lai():
    wit = list('ABCDEFG'); rows = []
    for u in json.load(open(os.path.join(P, 'historic/lai_ombre/units_v2.json'))):
        r = ['0'] * 7
        for k, g in enumerate(u['groups']):
            for c in g:
                if c in wit: r[wit.index(c)] = str(k)
        rows.append(r)
    return wit, rows


data = {k: load_num(v) for k, v in SRC.items()}
data["Lai de l'Ombre (Bedier 1913)"] = load_lai()
out = {}
for name, (wit, rows) in data.items():
    n = len(wit)
    def singular(r, i): return r[i] != '?' and sum(1 for k in range(n) if r[k] == r[i]) == 1
    p = [sum(singular(r, i) for r in rows) / max(1, sum(r[i] != '?' for r in rows)) for i in range(n)]
    pairs = []
    for i, j in itertools.combinations(range(n), 2):
        both = [r for r in rows if r[i] != '?' and r[j] != '?']
        if not both or p[i] * p[j] == 0: continue
        ret = sum(1 for r in both if singular(r, i) and singular(r, j)
                  and len({r[k] for k in range(n) if k not in (i, j) and r[k] != '?'}) == 1
                  and sum(1 for k in range(n) if k not in (i, j) and r[k] != '?') >= 2) / len(both)
        coi = sum(1 for r in both if r[i] == r[j] and sum(1 for k in range(n) if r[k] == r[i]) == 2
                  and sum(1 for k in range(n) if r[k] not in ('?', r[i])) >= 2) / len(both)
        pairs.append(dict(i=wit[i], j=wit[j], pp=p[i] * p[j], retention=ret / (p[i] * p[j]), coincidence=coi / (p[i] * p[j]),
                          retention_rate=ret, coincidence_rate=coi))
    pairs.sort(key=lambda d: -d['pp'])
    top = pairs[:5]
    k10 = pairs[:10]
    R10, C10 = sum(d['retention_rate'] for d in k10), sum(d['coincidence_rate'] for d in k10)
    Ra, Ca = sum(d['retention_rate'] for d in pairs), sum(d['coincidence_rate'] for d in pairs)
    multi = sum(1 for r in rows if len({v for v in r if v != '?'}) >= 3) / len(rows)
    out[name] = dict(n=n, pairs=len(pairs), retention_median=round(st.median(d['retention'] for d in pairs), 3),
                     coincidence_median=round(st.median(d['coincidence'] for d in pairs), 3),
                     units_with_3plus_readings=round(multi, 3),
                     top10_retention=round(R10, 5), top10_coincidence=round(C10, 5),
                     top10_C_over_R=round(C10 / R10, 3) if R10 else None,
                     all_retention=round(Ra, 4), all_coincidence=round(Ca, 4), all_C_over_R=round(Ca / Ra, 3) if Ra else None,
                     top5=[{k: (round(v, 4) if isinstance(v, float) else v) for k, v in d.items()} for d in top])
    print(name, {k: out[name][k] for k in ('pairs', 'units_with_3plus_readings', 'top10_retention', 'top10_coincidence', 'top10_C_over_R', 'all_retention', 'all_coincidence', 'all_C_over_R')})
    for d in top: print('   ', d['i'], d['j'], 'pp %.4f ret %.3f coi %.3f' % (d['pp'], d['retention'], d['coincidence']))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
