"""CHECK (2026-11-23): reruns of the frozen series-test rule (bedier_top.py via
bedier_sensitivity.classify: exact rooted MWIS min weight 2, B = 200 unit bootstrap seed 1, stable >= 0.90, contradicting
family weight >= 5) under alternative codings of the same extracted units. Analysis witnesses and the editor's tops are
held at the frozen run's (classification_a10). Not pre-registered; nothing frozen is modified.

Conditions
  printed      frozen coding (departures from the printed text); must reproduce classification_a10
  maj_all      B6: departures from the MAJORITY reading of each unit. Negative apparatus: every analysis witness not
               cited in a unit reads the lemma (all witnesses legible everywhere). If one cited group holds more than
               half of the legible witnesses, it becomes the reference: the lemma group (uncited + lemma sigla) and the
               other cited groups depart instead.
  maj_span     as maj_all, but a witness is legible in a unit only between the first and last leaf on which it is
               cited (proxy for fragmentary witnesses / lacunae; the apparatus does not record extent).
  nosp_strict  B2: printed coding with spelling-only departures removed (check_spelling.py, strict normalisation)
  nosp_broad   as above, broad (phonetic) normalisation
  nosp_union   broad, spelling-only if the printed lemma OR the edited-text line window says so (maximal removal)
For every condition: verdict per text, the size-matched random-division null of bedier_null.py (exact enumeration),
and the pooled statistic of the paper (sum of scores supported 2 / silent 1 / contradicted 0, exact convolution of the
13 per-text nulls, one-sided P(sum >= observed)); also the median over texts of the share of random divisions contradicted.
B10: on the printed run, the null restricted to divisions COMPATIBLE with the stable families (see compat_null).
Usage: check_bedier.py [--workers 12] [--conds printed,maj_all,...]
"""
import json, os, sys, argparse, collections, itertools, statistics, time
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bedier_sensitivity import classify, verdict, null_dist, crosses, SCORE, TEXTS, ROOT
import check_spelling as RS

OUT = os.path.join(ROOT, 'checks_20261123')
CONDS = ['printed', 'maj_all', 'maj_span', 'nosp_strict', 'nosp_broad', 'nosp_union']
SETTINGS = {s['text_id']: s for s in json.load(open(os.path.join(ROOT, 'settings.json')))}


def base(t):
    d = json.load(open(os.path.join(ROOT, 'classification_a10', t + '.json')))
    U = json.load(open(os.path.join(ROOT, 'runs', t, 'units.json')))
    W = list(d['analysis_witnesses']); Ws = frozenset(W)
    return d, U, W, Ws, [frozenset(b) for b in d['editor_top']]


def fam_printed(U, Ws, drop=None):
    F = []
    for i, u in enumerate(U):
        gs = u['groups'][1:]
        fs = [frozenset(x for x in g if x in Ws) for j, g in enumerate(gs) if not (drop and drop[i][j])]
        fs = [g for g in fs if len(g) >= 2 and g != Ws]
        if fs: F.append(fs)
    return F


def fam_majority(U, W, Ws, span):
    rng = {}
    if span:
        for u in U:
            for g in u['groups']:
                for x in g:
                    if x in Ws:
                        lo, hi = rng.get(x, (u['leaf'], u['leaf'])); rng[x] = (min(lo, u['leaf']), max(hi, u['leaf']))
    F, flips, units = [], 0, 0
    for u in U:
        leg = frozenset(w for w in W if not span or (w in rng and rng[w][0] <= u['leaf'] <= rng[w][1]))
        cited = [frozenset(x for x in g if x in leg) for g in u['groups'][1:]]
        cited = [g for g in cited if g]
        lemma = leg - frozenset().union(*cited) if cited else leg
        big = [g for g in cited if 2 * len(g) > len(leg)]
        if big:
            flips += 1; ref = big[0]
            dep = [lemma] + [g for g in cited if g is not ref]
        else:
            dep = cited
        fs = [g for g in dep if len(g) >= 2 and g != Ws and g != leg]
        if fs: F.append(fs); units += 1
    return F, flips


def spelling_mask(t, U, broad, union):
    L = RS.Leaves(t, SETTINGS[t])
    return [RS.classify_unit(u, L, broad, union) for u in U]


def job(args):
    cond, t = args; t0 = time.time()
    d, U, W, Ws, top = base(t)
    extra = {}
    if cond == 'printed': F = fam_printed(U, Ws)
    elif cond.startswith('maj'):
        F, flips = fam_majority(U, W, Ws, span=cond == 'maj_span'); extra['units_flipped'] = flips
    else:
        mask = spelling_mask(t, U, broad=cond != 'nosp_strict', union=cond == 'nosp_union')
        n_dep = sum(1 for i, u in enumerate(U) for j, g in enumerate(u['groups'][1:]) if len(frozenset(g) & Ws) >= 2)
        n_sp = sum(1 for i, u in enumerate(U) for j, g in enumerate(u['groups'][1:]) if len(frozenset(g) & Ws) >= 2 and mask[i][j])
        extra.update(shared_departures=n_dep, spelling_only_removed=n_sp)
        F = fam_printed(U, Ws, drop=mask)
    v, stable, bad = classify(F, top)
    br = []
    for b in top:
        if len(b) < 2: br.append(dict(branch=''.join(sorted(b)), status='single_witness')); continue
        S = {f for f, _ in stable}
        if b in S: br.append(dict(branch=''.join(sorted(b)), status='supported'))
        else:
            cr = [(''.join(sorted(f)), w) for f, w in stable if w >= 5 and crosses(f, b)]
            br.append(dict(branch=''.join(sorted(b)), status='contradicted' if cr else 'silent', by=cr or None))
    return dict(cond=cond, text=t, verdict=v, frozen=d['classification'], units=len(F), frozen_units=d['units_with_shared_departure'],
                stable=[[''.join(sorted(f, key=d['stemma_witnesses'].index)), w] for f, w in sorted(stable, key=lambda x: -x[1])],
                branches=br, null=null_dist(W, top, stable), not_optimal=bad, secs=round(time.time() - t0), **extra)


def pooled(rows):
    used = [r for r in rows if r['null']]
    out = {}
    for name, fun in [('score_sum', lambda v: SCORE[v]), ('n_supported', lambda v: v == 'supported'),
                      ('n_not_contradicted', lambda v: v != 'contradicted')]:
        D = {0: 1.0}
        for r in used:
            nd = collections.defaultdict(float)
            for s0, p0 in D.items():
                for v, q in r['null'].items(): nd[s0 + fun(v)] += p0 * q
            D = nd
        obs = sum(fun(r['verdict']) for r in used)
        out[name] = dict(observed=obs, null_mean=sum(k * v for k, v in D.items()), p_one_sided=sum(v for k, v in D.items() if k >= obs - 1e-9))
    out['median_share_random_contradicted'] = statistics.median(r['null'].get('contradicted', 0) for r in used)
    out['n_texts'] = len(used)
    out['editor_counts'] = dict(collections.Counter(r['verdict'] for r in rows))
    return out


def compat_null(t, stable, which):
    """B10. Divisions with the editor's branch sizes (drawn from the analysis witnesses, as in bedier_null.py) that cross
    no 'lower' stable family. which='all': every stable family of the frozen optimum (then 'contradicted' is impossible
    by construction); which='nonmax': only stable families strictly inside another stable family (the maximal ones,
    which are themselves candidate tops, are not imposed)."""
    d, U, W, Ws, top = base(t)
    if len(top) != 2: return None
    fams = [f for f, _ in stable]
    if which == 'nonmax': cons = [f for f in fams if any(f < g for g in fams)]
    else: cons = fams
    s1, s2 = sorted([len(top[0]), len(top[1])]); dist = collections.Counter(); seen = set(); n_all = 0
    for c1 in itertools.combinations(W, s1):
        rest = [w for w in W if w not in c1]
        for c2 in itertools.combinations(rest, s2):
            k = frozenset([frozenset(c1), frozenset(c2)])
            if k in seen: continue
            seen.add(k); n_all += 1
            b1, b2 = frozenset(c1), frozenset(c2)
            if any(crosses(f, b) for f in cons for b in (b1, b2)): continue
            dist[verdict([b1, b2], stable)] += 1
    N = sum(dist.values())
    ed = verdict(top, stable)
    ed_compat = not any(crosses(f, b) for f in cons for b in top)
    e = SCORE[ed]
    res = dict(n_divisions=n_all, n_compatible=N, constraints=[''.join(sorted(f)) for f in cons], editor=ed,
               editor_compatible=ed_compat, counts=dict(dist))
    if N:
        res['shares'] = {k: v / N for k, v in dist.items()}
        res['n_better'] = sum(v for k, v in dist.items() if SCORE[k] > e)
        res['n_tied'] = sum(v for k, v in dist.items() if SCORE[k] == e)
        res['n_worse'] = sum(v for k, v in dist.items() if SCORE[k] < e)
        res['p_one_sided'] = (res['n_better'] + res['n_tied']) / N
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--conds', default=','.join(CONDS)); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    conds = a.conds.split(',')
    jobs = sorted([(c, t) for c in conds for t in TEXTS], key=lambda x: x[1] in ('claudianus_maior', 'ennodius_opera', 'historia_brittonum'), reverse=True)
    with Pool(a.workers) as p:
        rows = []
        for r in p.imap_unordered(job, jobs):
            rows.append(r); print(f"{r['cond']:12s} {r['text']:32s} {r['verdict']:14s} (frozen {r['frozen']}) units {r['units']} {r['secs']}s", flush=True)
    res = {}
    for c in conds:
        rr = sorted([r for r in rows if r['cond'] == c], key=lambda r: TEXTS.index(r['text']))
        res[c] = dict(texts=rr, pooled=pooled(rr))
        json.dump(res[c], open(os.path.join(OUT, f'cond_{c}.json'), 'w'), indent=1)
    if 'printed' in res:
        comp = {}
        for r in res['printed']['texts']:
            st = [(frozenset(k), w) for k, w in r['stable']]
            comp[r['text']] = {w: compat_null(r['text'], st, w) for w in ('all', 'nonmax')}
        json.dump(comp, open(os.path.join(OUT, 'compat_null.json'), 'w'), indent=1)
    print(json.dumps({c: res[c]['pooled'] for c in res}, indent=1))


if __name__ == '__main__':
    main()
