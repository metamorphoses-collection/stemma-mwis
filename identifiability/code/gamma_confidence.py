r"""A lower confidence bound for the margin gamma that remains valid when the stemma is chosen from the data.

Theorem gammaconf: let w(s) ~ Binomial(L, m_s) for every non-trivial split s (units i.i.d.; a unit may propose several
splits). With one-sided Clopper-Pearson bounds at level delta/(2N), one lower and one upper per split, N = 2^(n-1) - n - 1, lower bounds l_s and upper bounds
u_s hold for every split simultaneously with probability >= 1 - delta (union bound). On that event, for EVERY stemma T,
  gamma(T) >= gamma_lo(T) := min over compatible S != T of (l(T\S) - u(S\T)) / |T\S|,
computed by the bisection of Lemma gamma with weights l on T and u elsewhere. Since the event does not depend on T,
the bound holds for the estimated stemma T_hat as well. For large n, splits never observed share the upper bound u_0
and are replaced by the correction - u_0 (1 + n - 3 - |T|) (Remark in the paper); here n <= 10, so all splits are
enumerated and no correction is needed.

Parts: (1) coverage and certification rate on the ten-witness regimes with exact m (multinomial sampling: each unit
proposes at most one split there); (2) the Lai, departures from the printed text, T_hat = {AB, CG, DE, DEF}.
Usage: gamma_confidence.py OUT.json"""
import sys, os, json, itertools, math
import numpy as np
from collections import Counter
from scipy.stats import beta
from ortools.sat.python import cp_model

DELTA = 0.05
rng = np.random.default_rng(20261118)


def all_splits(n):
    X = frozenset(range(n)); out = []
    for r in range(2, n - 1):
        for c in itertools.combinations(range(1, n), r):     # side not containing witness 0
            out.append(frozenset(c))
    return X, out


def conflict(A, B, X):
    return bool(A & B and A - B and B - A and X - (A | B))


def cp_bounds(k, L, d):
    lo = 0.0 if k == 0 else float(beta.ppf(d, k, L - k + 1))
    hi = 1.0 if k == L else float(beta.ppf(1 - d, k + 1, L - k))
    return lo, hi


class Solver:
    """Max-weight compatible set. Splits are keyed by their side without witness 0, so compatibility is laminarity
    of subsets of {1..n-1}; solved exactly by the O(3^n) dynamic programme over laminar families (Theorem FPT).
    best_cpsat is the MWIS formulation, kept as a cross-check."""
    def __init__(self, X, splits):
        self.V = splits; self.X = X; self.n = len(X)
        self.E = None

    def best(self, wt):
        k = self.n - 1; full = (1 << k) - 1
        bit = lambda s: sum(1 << (i - 1) for i in s)
        W = {}
        for v in self.V:
            W[bit(v)] = max(0.0, wt[v])
        H = [0.0] * (1 << k); P = [0.0] * (1 << k); Hc = [None] * (1 << k); Pc = [None] * (1 << k)
        for S in range(1, 1 << k):
            if S & (S - 1) == 0:
                H[S] = P[S] = 0.0; Pc[S] = (S, 0); continue
            low = S & -S; rest = S ^ low
            bestF, arg = -1.0, None
            T = rest
            while True:                                   # proper subsets T' = low | T, T' != S
                Tp = low | T
                if Tp != S:
                    val = H[Tp] + P[S ^ Tp]
                    if val > bestF: bestF, arg = val, Tp
                if T == 0: break
                T = (T - 1) & rest
            add = W.get(S, 0.0) if S != full else 0.0
            H[S] = bestF + add; Hc[S] = (arg, add > 0)
            if H[S] >= bestF: P[S], Pc[S] = H[S], (S, 0)
            else: P[S], Pc[S] = bestF, (arg, 1)
        chosen = set()
        def walkH(S):
            if S & (S - 1) == 0: return
            arg, take = Hc[S]
            if take: chosen.add(S)
            walkH(arg); walkP(S ^ arg)
        def walkP(S):
            if S == 0: return
            blk, _ = Pc[S]
            if blk == S: walkH(S)
            else: walkH(blk); walkP(S ^ blk)
        walkH(full)
        inv = {bit(v): v for v in self.V}
        sel = {inv[c] for c in chosen if c in inv}
        return sum(max(0.0, wt[v]) for v in sel), sel

    def best_cpsat(self, wt):
        if self.E is None:
            self.E = [(i, j) for i, j in itertools.combinations(range(len(self.V)), 2) if conflict(self.V[i], self.V[j], self.X)]
        m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in self.V]
        for i, j in self.E: m.Add(x[i] + x[j] <= 1)
        m.Maximize(sum(int(round(1e9 * wt[v])) * x[k] for k, v in enumerate(self.V)))
        s = cp_model.CpSolver(); s.parameters.num_workers = 8; s.Solve(m)
        return s.ObjectiveValue() / 1e9, {v for k, v in enumerate(self.V) if s.Value(x[k])}


def gamma_of(solver, T, wT, wO, tol=1e-7):
    """gamma for weights wT on T and wO elsewhere, by bisection (Lemma gamma)."""
    def T_is_max(g):
        wt = {v: (wT[v] - g if v in T else wO[v]) for v in solver.V}
        opt, _ = solver.best(wt)
        return opt <= sum(wt[v] for v in T) + 1e-9
    lo, hi = -1.0, 1.0
    if not T_is_max(lo): return -1.0
    while hi - lo > tol:
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if T_is_max(mid) else (lo, mid)
    return lo


def T_is_max_at(solver, T, wT, wO, g):
    """Lemma gamma: T is a maximum of the shifted weights iff g <= gamma(wT, wO); one MWIS."""
    wt = {v: (wT[v] - g if v in T else wO[v]) for v in solver.V}
    opt, _ = solver.best(wt)
    return opt <= sum(wt[v] for v in T) + 1e-9


def bounds(solver, cnt, L, n):
    d = DELTA / (2 * (2 ** (n - 1) - n - 1))
    b = {v: cp_bounds(cnt.get(v, 0), L, d) for v in solver.V}
    return {v: b[v][0] for v in solver.V}, {v: b[v][1] for v in solver.V}


def certify(solver, cnt, L, T, n):
    N = 2 ** (n - 1) - n - 1; d = DELTA / (2 * N)
    b = {v: cp_bounds(cnt.get(v, 0), L, d) for v in solver.V}
    g_lo = gamma_of(solver, T, {v: b[v][0] for v in solver.V}, {v: b[v][1] for v in solver.V})
    g_hat = gamma_of(solver, T, {v: cnt.get(v, 0) / L for v in solver.V}, {v: cnt.get(v, 0) / L for v in solver.V})
    return g_lo, g_hat, b


res = {'delta': DELTA}
# ---- check: DP equals CP-SAT (n = 7, random weights; n = 10, two instances) ----
import random as _r; _r.seed(5); dp_ok = True
for nn, reps in ((7, 30), (10, 2)):
    Xc, Sc = all_splits(nn); sc = Solver(Xc, Sc)
    for _ in range(reps):
        wt = {v: (_r.random() if _r.random() < 0.3 else 0.0) for v in Sc}
        a, sel = sc.best(wt); b = sc.best_cpsat(wt)[0]
        laminar = all(not conflict(x, y, Xc) for x, y in itertools.combinations(sel, 2))
        dp_ok &= abs(a - b) < 1e-6 and laminar
res['dp_equals_cpsat'] = bool(dp_ok); print('DP = CP-SAT:', dp_ok, flush=True)
assert dp_ok
# ---- part 1: ten witnesses, exact m ----
X, S10 = all_splits(10); sol10 = Solver(X, S10)
key = lambda s: frozenset(int(t) for t in s.split('|'))
for reg in json.load(open('data/run2/exact_population.json')):
    m = {key(k): v for k, v in eval(reg['m']).items()} if isinstance(reg['m'], str) else {key(k): v for k, v in reg['m'].items()}
    m = {(s if 0 not in s else X - s): v for s, v in m.items()}
    T = {key(s) for s in (eval(reg['population_optimum']) if isinstance(reg['population_optimum'], str) else reg['population_optimum'])} if str(reg['population_optimum_is_true']) == 'True' else None
    if T is None: continue                     # long-branch: the true set is not the maximum, gamma < 0
    T = {(s if 0 not in s else X - s) for s in T}
    g_true = float(reg['gamma']); g_check = gamma_of(sol10, T, m, m)
    assert abs(g_check - g_true) < 1e-5, (g_check, g_true)
    p = np.array([m.get(v, 0.0) for v in S10]); p0 = 1 - p.sum()
    rows = []
    for L in (1600, 6400, 25600, 102400, 409600):
        cover = cert = cert_true = 0; R = 40
        for _ in range(R):
            c = rng.multinomial(L, np.append(p, p0))[:-1]
            cnt = {v: int(k) for v, k in zip(S10, c) if k}
            lo, hi = bounds(sol10, cnt, L, 10)
            cover += not T_is_max_at(sol10, T, lo, hi, g_true + 1e-7)          # gamma_lo(T) <= gamma(T)
            _, That = sol10.best({v: cnt.get(v, 0) / L for v in S10})
            ok = len(That) == 7 and T_is_max_at(sol10, That, lo, hi, 1e-7)       # binary T_hat only; gamma_lo > 0
            cert += ok; cert_true += ok and That == T
        rows.append(dict(L=L, reps=R, coverage=cover / R, certified=cert / R, certified_and_true=cert_true / R))
        print(reg['name'], rows[-1], flush=True)
    res[reg['name']] = dict(gamma=g_true, rows=rows)

# ---- part 2: the Lai ----
W = 'ABCDEFG'; ix = {c: k for k, c in enumerate(W)}
X7, S7 = all_splits(7); sol7 = Solver(X7, S7)
U = json.load(open(os.path.expanduser('~/MC/projects/priscus_excerpta/data/historic/lai_ombre/units.json')))
# hybrid coding: departure classes (non-lemma readings) of sizes 2..5, each identified with its unordered split;
# a unit counts each split at most once (two departure classes that are complements give one split)
cnt = Counter()
for u in U:
    props = set()
    for k, g in enumerate(u['groups']):
        s = frozenset(ix[c] for c in g if c in W)
        if k != u['lemma_group'] and 2 <= len(s) <= 5: props.add(s if 0 not in s else X7 - s)
    for s in props: cnt[s] += 1
L = len(U)
T = {frozenset(ix[c] for c in t) for t in ('CG', 'DE', 'DEF')} | {X7 - frozenset({0, 1})}
_, That = sol7.best({v: cnt.get(v, 0) / L for v in S7})
assert len(T) == 4 and That == T                 # binary on 7 witnesses
g_lo, g_hat, b = certify(sol7, cnt, L, T, 7)
name = lambda s: ''.join(W[i] for i in sorted(s))
res['lai'] = dict(L=L, T=[name(s) for s in T], T_is_optimum=That == T, gamma_hat=g_hat, gamma_lo=g_lo,
                  supports={name(s): cnt.get(s, 0) for s in T},
                  bounds={name(s): [round(b[s][0], 5), round(b[s][1], 5)] for s in T},
                  units_needed_per_split=math.log(2 * 56 / DELTA))
print('Lai', res['lai'])
json.dump(res, open(sys.argv[1], 'w'), indent=1, default=str)

# ---- part 3: how much more text the Lai would need, at the same proportions ----
def lai_lo(k):
    c2 = {v: int(round(k * x)) for v, x in cnt.items()}
    return certify(sol7, c2, int(round(k * L)), T, 7)[0]
ks = [1, 1.5, 2, 3, 4, 6, 8]
res['lai']['scaled'] = {k: lai_lo(k) for k in ks}
print('Lai scaled', res['lai']['scaled'])
json.dump(res, open(sys.argv[1], 'w'), indent=1, default=str)

# ---- check C22: fails (non-zero exit) if any expected property fails ----
fails = []
if not res['dp_equals_cpsat']: fails.append('DP != CP-SAT')
for reg in ('clean', 'noisy'):
    for r in res[reg]['rows']:
        if r['coverage'] < 1 - DELTA: fails.append(f'{reg} L={r["L"]}: coverage {r["coverage"]} below 1 - delta')
        if r['certified'] != r['certified_and_true']: fails.append(f'{reg} L={r["L"]}: certified a wrong stemma')
sc = res['lai']['scaled']
if not all(sc[a] <= sc[b] for a, b in zip(list(sc), list(sc)[1:])): fails.append('Lai bound not monotone in length')
print('C22', 'ALL PASS' if not fails else fails); sys.exit(1 if fails else 0)
