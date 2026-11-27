"""EXPLORATORY (not pre-registered; written 2026-11-22): sensitivity of the series test
(bedier_top.py, frozen rule) to extraction error.

For each of the 14 primary texts the extracted families (groups of >= 2 analysis witnesses sharing a departure, after
the frozen witness filter) are perturbed, and the frozen verdict is recomputed in full: exact rooted MWIS (min weight 2),
B = 200 unit bootstrap (seed 1), stable >= 0.90, contradicting family weight >= 5. The analysis witnesses are held at
those of the frozen run; the 200-unit inclusion gate is recorded but not applied.

Scenarios (N = number of extracted family occurrences of the text):
  C0   control: no perturbation, bootstrap seed varied (Monte Carlo noise of the stability estimate itself)
  D16  random misses: each family occurrence deleted with p = 0.16 (measured recall 0.84)
  D30  random misses, p = 0.30 (pessimistic)
  A10  spurious groups, false share 0.10 of the result (measured precision >= 0.90): half are groups of a
       neighbouring page (leaf +-1) copied into a unit, half are misread sigla (one siglum substituted, added or
       dropped, replacing the true group in place)
  A20  as A10, false share 0.20
  M    D16 then A10 (the measured profile combined)
  BX   biased misses, 16% overall, groups that cross an editor branch 4x likelier to be missed than others
  BA   biased misses, 16% overall, groups inside an editor branch 4x likelier to be missed
  RX   missed rival readings restored: 0.16/0.84 N extra occurrences, all copied from groups that cross an editor
       branch (as if every missed group had been a rival reading) -- the worst case for 'supported' and 'silent'
  RA   missed agreeing readings restored: 0.16/0.84 N extra occurrences, all copied from groups inside an editor
       branch -- the worst case for 'silent' and 'contradicted'
Each (scenario, text, replicate) also stores the exact size-matched random-top null of bedier_null.py on that
replicate's stable families, so the pooled 'editors vs random tops' p-value can be recomputed per replicate.

Usage: bedier_sensitivity.py --out DIR [--reps 200] [--workers 16] [--scenarios C0,D16,...]
       bedier_sensitivity.py --baseline --out DIR      (reproduce the frozen classification first)
"""
import json, os, sys, argparse, random, collections, itertools, time, zlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rooted_mwis import mwis, weigh

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'bedier')
TEXTS = ['isidorus_chronica_maiora', 'claudianus_raptus_proserpinae', 'prosper_epitoma_chronicon', 'historia_brittonum',
         'isidorus_historia_gothorum', 'fortunatus_vita_martini', 'claudianus_maior', 'sidonius_carmina',
         'ennodius_opera', 'paulus_historia_romana', 'victor_vitensis_historia', 'jordanes_romana',
         'eugippius_vita_severini', 'marcellinus_comes_chronicon']
SCEN = ['C0', 'D16', 'D30', 'A10', 'A20', 'M', 'BX', 'BA', 'RX', 'RA']
STABLE, CROSS_W, MINW, B, MIN_UNITS = 0.90, 5, 2, 200, 200
SCORE = {'supported': 2, 'partly_silent': 1, 'contradicted': 0}


def crosses(a, b): return bool(a & b) and not (a <= b or b <= a)


def load(t):
    d = json.load(open(os.path.join(ROOT, 'classification_a10', t + '.json')))
    U = json.load(open(os.path.join(ROOT, 'runs', t, 'units.json')))
    W = list(d['analysis_witnesses']); Ws = frozenset(W)
    F, L = [], []
    for u in U:
        fs = [frozenset(x for x in g if x in Ws) for g in u['groups'][1:]]
        fs = [g for g in fs if len(g) >= 2 and g != Ws]
        if fs: F.append(fs); L.append(u['leaf'])
    top = [frozenset(b) for b in d['editor_top']]
    return dict(d=d, W=W, Ws=Ws, F=F, L=L, top=top, order=d['stemma_witnesses'])


def verdict(top, stable):
    S = {f for f, _ in stable}; st = []
    for b in top:
        if len(b) < 2: continue
        if b in S: st.append('supported')
        elif any(w >= CROSS_W and crosses(f, b) for f, w in stable): st.append('contradicted')
        else: st.append('silent')
    if not st: return None
    return 'contradicted' if 'contradicted' in st else ('supported' if all(s == 'supported' for s in st) else 'partly_silent')


def classify(F, top, seed=1):
    Wt = weigh(F, MINW); sel, obj, status = mwis(Wt)
    rng = random.Random(seed); stab = collections.Counter(); bad = status != 'OPTIMAL'
    for _ in range(B):
        Fb = [F[rng.randrange(len(F))] for _ in F]
        sb, _, s2 = mwis(weigh(Fb, MINW)); stab.update(sb); bad |= s2 != 'OPTIMAL'
    stable = [(f, Wt[f]) for f in sel if stab[f] / B >= STABLE]
    return verdict(top, stable), stable, bad


def null_dist(W, top, stable):
    if len(top) != 2: return None
    s1, s2 = sorted([len(top[0]), len(top[1])]); dist = collections.Counter(); seen = set()
    for c1 in itertools.combinations(W, s1):
        rest = [w for w in W if w not in c1]
        for c2 in itertools.combinations(rest, s2):
            k = frozenset([frozenset(c1), frozenset(c2)])
            if k in seen: continue
            seen.add(k); dist[verdict([frozenset(c1), frozenset(c2)], stable)] += 1
    N = sum(dist.values()); return {k: v / N for k, v in dist.items()}


def kind(g, top):
    br = [b for b in top if len(b) >= 2]
    if any(crosses(g, b) for b in br): return 'X'
    if any(g <= b for b in br): return 'A'
    return 'O'


def delete(F, rng, p_of):
    out = []
    for fs in F:
        k = [g for g in fs if rng.random() >= p_of(g)]
        if k: out.append(k)
    return out


def biased_p(F, top, target, ratio, cls):
    n = sum(len(fs) for fs in F); nb = sum(1 for fs in F for g in fs if kind(g, top) == cls)
    if nb == 0: return lambda g: target
    c = target * n / (ratio * nb + (n - nb))  # p = c for others, ratio*c for the favoured class, capped at 1
    if ratio * c > 1: c = (target * n - nb) / (n - nb); return lambda g: 1.0 if kind(g, top) == cls else c
    return lambda g: ratio * c if kind(g, top) == cls else c


def corrupt(g, W, rng):
    for _ in range(50):
        h = set(g); op = rng.choice('sad')
        if op == 's':
            h.discard(rng.choice(sorted(h))); h.add(rng.choice([w for w in W if w not in g] or W))
        elif op == 'a':
            cand = [w for w in W if w not in g]
            if not cand: continue
            h.add(rng.choice(cand))
        elif len(h) > 2: h.discard(rng.choice(sorted(h)))
        else: continue
        h = frozenset(h)
        if len(h) >= 2 and h != frozenset(W) and h != g: return h
    return None


def add_spurious(F, L, W, q, rng):
    F = [list(fs) for fs in F]; n = sum(len(fs) for fs in F); n_false = round(q / (1 - q) * n)
    byleaf = collections.defaultdict(list)
    for i, l in enumerate(L): byleaf[l].append(i)
    done = 0; tries = 0
    while done < n_false and tries < 20 * n_false + 100:
        tries += 1; i = rng.randrange(len(F))
        if rng.random() < 0.5:  # a correct group of the neighbouring page lands in this unit
            nb = byleaf.get(L[i] + rng.choice((-1, 1)))
            if not nb: continue
            g = rng.choice(F[rng.choice(nb)]); F[i].append(g); done += 1
        else:  # a misread siglum: the true group is replaced by a corrupted one
            j = rng.randrange(len(F[i])); h = corrupt(F[i][j], W, rng)
            if h is None: continue
            F[i][j] = h; done += 1
    return F


def restore(F, top, cls, rng):
    pool = [g for fs in F for g in fs if kind(g, top) == cls]
    F = [list(fs) for fs in F]
    if not pool: return F, 0
    n_add = round(0.16 / 0.84 * sum(len(fs) for fs in F))
    for _ in range(n_add): F[rng.randrange(len(F))].append(rng.choice(pool))
    return F, n_add


_CACHE = {}


def task(args):
    sc, t, rep = args
    if t not in _CACHE: _CACHE[t] = load(t)
    D = _CACHE[t]; F, top, W = D['F'], D['top'], D['W']
    rng = random.Random(zlib.crc32(f"{sc}|{t}|{rep}".encode())); seed = 1; extra = {}
    if sc == 'C0': seed = 1000 + rep
    elif sc == 'D16': F = delete(F, rng, lambda g: 0.16)
    elif sc == 'D30': F = delete(F, rng, lambda g: 0.30)
    elif sc == 'A10': F = add_spurious(F, D['L'], W, 0.10, rng)
    elif sc == 'A20': F = add_spurious(F, D['L'], W, 0.20, rng)
    elif sc == 'M':
        keep = [(fs2, l) for fs, l in zip(F, D['L']) for fs2 in [[g for g in fs if rng.random() >= 0.16]] if fs2]
        F = add_spurious([k for k, _ in keep], [l for _, l in keep], W, 0.10, rng)
    elif sc == 'BX': F = delete(F, rng, biased_p(F, top, 0.16, 4, 'X'))
    elif sc == 'BA': F = delete(F, rng, biased_p(F, top, 0.16, 4, 'A'))
    elif sc == 'RX': F, extra['n_add'] = restore(F, top, 'X', rng)
    elif sc == 'RA': F, extra['n_add'] = restore(F, top, 'A', rng)
    t0 = time.time(); v, stable, bad = classify(F, top, seed)
    key = lambda f: ''.join(sorted(f, key=D['order'].index))
    return dict(scenario=sc, text=t, rep=rep, verdict=v, n_units=len(F), gate_200=len(F) >= MIN_UNITS,
                non_optimal=bad, stable=[[key(f), w] for f, w in stable], null=null_dist(W, top, stable),
                secs=round(time.time() - t0, 2), **extra)


def baseline(out):
    rows = []
    for t in TEXTS:
        D = load(t); v, stable, bad = classify(D['F'], D['top'], 1)
        key = lambda f: ''.join(sorted(f, key=D['order'].index))
        frozen = D['d']; fs = sorted([key(f), w] for f, w in stable); ff = sorted([k, w] for k, w, s in frozen['stable_families'])
        rows.append(dict(text=t, verdict=v, frozen=frozen['classification'], match=(v == frozen['classification'] and fs == ff),
                         units=len(D['F']), frozen_units=frozen['units_with_shared_departure'], stable=fs,
                         branches=frozen['branches'], null=null_dist(D['W'], D['top'], stable)))
        print(t, v, frozen['classification'], rows[-1]['match'], len(D['F']), frozen['units_with_shared_departure'])
    json.dump(rows, open(os.path.join(out, 'baseline.json'), 'w'), indent=1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--reps', type=int, default=200)
    ap.add_argument('--workers', type=int, default=16); ap.add_argument('--scenarios', default=','.join(SCEN))
    ap.add_argument('--baseline', action='store_true'); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.baseline: baseline(a.out); return
    path = os.path.join(a.out, 'replicates.jsonl'); done = set()
    if os.path.exists(path):
        for line in open(path): r = json.loads(line); done.add((r['scenario'], r['text'], r['rep']))
    jobs = [(sc, t, r) for t in TEXTS for sc in a.scenarios.split(',') for r in range(a.reps) if (sc, t, r) not in done]
    print(len(jobs), 'jobs', flush=True)
    import multiprocessing as mp
    t0 = time.time()
    with mp.Pool(a.workers) as pool, open(path, 'a') as fh:
        for i, r in enumerate(pool.imap_unordered(task, jobs, chunksize=4)):
            fh.write(json.dumps(r) + '\n')
            if i % 500 == 0: fh.flush(); print(i, round(time.time() - t0), 's', flush=True)
    print('done', round(time.time() - t0), 's')


if __name__ == '__main__':
    main()
