"""Exact expected supports m for the ten-witness regimes of sim_identifiability.py (all 2^18 innovation patterns of
the irreversible model, mixed with the contamination component exactly as simulated: with probability eta the unit's
reading set is a true clade drawn uniformly from the 8 non-trivial clades with one witness, drawn uniformly from the
10, toggled in or out). From m: exact global/local margins, the crossing pair that violates local separation, the
margin gamma (bisection, one exact MWIS per step, bracket [-1, 1], 40 steps), the heaviest competitor when gamma < 0,
and both bounds at the simulated lengths. Usage: exact_population.py RUN2_JSON OUT_JSON"""
import sys, json, math, itertools
import numpy as np
from ortools.sat.python import cp_model

sys.path.insert(0, __file__.rsplit('/', 1)[0])
n = 10
clades = [frozenset({i}) for i in range(n)] + [frozenset(c) for c in
          [{0, 1}, {0, 1, 2}, {3, 4}, {0, 1, 2, 3, 4}, {5, 6}, {5, 6, 7}, {8, 9}, {5, 6, 7, 8, 9}]]
X = frozenset(range(n))
canon = lambda Y: frozenset(Y) if 0 not in Y else X - frozenset(Y)
nontriv = lambda Y: 2 <= len(Y) <= n - 2
def compatible(A, B):
    Ac, Bc = X - A, X - B
    return not (A & B and A & Bc and Ac & B and Ac & Bc)
true = {canon(c) for c in clades if nontriv(c)}
masks = np.array([sum(1 << i for i in c) for c in clades], dtype=np.int64)


def exact_m(rates, eta):
    E = len(clades); r = np.array(rates)
    pats = np.arange(1 << E, dtype=np.int64)
    bits = ((pats[:, None] >> np.arange(E)) & 1).astype(bool)
    prob = np.prod(np.where(bits, r, 1 - r), axis=1)
    Y = np.zeros(len(pats), dtype=np.int64)
    for k in range(E): Y |= np.where(bits[:, k], masks[k], 0)
    m = {}
    for y, pr in zip(*np.unique(Y, return_inverse=False, return_counts=False), ) if False else []: pass
    acc = {}
    np.add.at(acc.setdefault('x', np.zeros(1 << n)), Y, prob)
    tot = acc['x']
    for y in np.flatnonzero(tot):
        S = frozenset(i for i in range(n) if y >> i & 1)
        if nontriv(S): s = canon(S); m[s] = m.get(s, 0) + (1 - eta) * tot[y]
    tc = [c for c in clades if nontriv(c)]
    for c in tc:
        for w in range(n):
            S = set(c) ^ {w}
            if nontriv(S): s = canon(S); m[s] = m.get(s, 0) + eta / (len(tc) * n)
    return m


def mwis_value(wts):
    V = [v for v, x in wts.items() if x > 0]
    md = cp_model.CpModel(); x = {v: md.NewBoolVar('') for v in V}
    for a, b in itertools.combinations(V, 2):
        if not compatible(a, b): md.Add(x[a] + x[b] <= 1)
    SC = 10 ** 9
    md.Maximize(sum(int(round(wts[v] * SC)) * x[v] for v in V))
    s = cp_model.CpSolver(); s.parameters.num_workers = 8; s.Solve(md)
    return s.ObjectiveValue() / SC, frozenset(v for v in V if s.Value(x[v]))


def gamma(m):
    def ok(g):
        sh = {v: x - (g if v in true else 0) for v, x in m.items()}
        val, _ = mwis_value(sh)
        return val <= sum(sh.get(v, 0) for v in true) + 1e-9
    lo, hi = -1.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2; lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


run = json.load(open(sys.argv[1])); out = []
for reg in run['regimes']:
    m = exact_m(reg['rates'], reg['eta'])
    b = min(m.get(s, 0) for s in true); a = max(v for s, v in m.items() if s not in true)
    pairs = [(m[s], m.get(t, 0), s, t) for s in m if s not in true for t in true if not compatible(s, t)]
    worst = max(pairs, key=lambda p: p[0] - p[1])
    g = gamma(m)
    val, S = mwis_value(m)
    rows = []
    for row in reg['rows']:
        L = row['L']
        pw = sum((1 - (math.sqrt(mt) - math.sqrt(ms)) ** 2) ** L for ms, mt, _, _ in pairs if mt > ms) \
            + sum(1 for ms, mt, _, _ in pairs if mt <= ms) + sum((1 - m.get(t, 0)) ** L for t in true)
        bern = sum(2 * math.exp(-L * g * g / (8 * v + 4 * g / 3)) for v in m.values()) if g > 0 else None
        rows.append(dict(L=L, recovery=row['recovery'], pairwise_bound=min(1.0, pw) if all(mt > ms for ms, mt, _, _ in pairs) else None,
                         bernstein_bound=min(1.0, bern) if bern is not None else None))
    fmt = lambda s: '|'.join(map(str, sorted(s)))
    out.append(dict(name=reg['name'], K=sum(1 for s in m if s not in true), min_true=b, max_false=a,
                    local_margin=worst[1] - worst[0], worst_pair=dict(false=fmt(worst[2]), m_false=worst[0], true=fmt(worst[3]), m_true=worst[1]),
                    gamma=g, population_optimum_is_true=S == true, population_optimum_value=val, true_value=sum(m.get(t, 0) for t in true),
                    population_optimum=[fmt(s) for s in S], rows=rows, m={fmt(s): v for s, v in sorted(m.items(), key=lambda kv: -kv[1])}))
    print(reg['name'], 'K', out[-1]['K'], 'min_true %.5f max_false %.5f local %.5f gamma %.5f opt_true %s' % (b, a, out[-1]['local_margin'], g, S == true),
          out[-1]['worst_pair'], [(r['L'], r['pairwise_bound'], r['bernstein_bound']) for r in rows])
json.dump(out, open(sys.argv[2], 'w'), indent=1)
