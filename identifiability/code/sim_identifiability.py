"""Simulation checks for the separation lemma and the finite-sample bound (paper_stemma_identifiability).

Model: rooted binary stemma, binary units, irreversible innovation (reading 0 -> 1 on edge e with prob pi_e,
no reversion), optional contamination (with prob eta a unit's reading set is a true clade with one witness
moved in or out). Each unit proposes the split Y|X\\Y of its 1-readers if non-trivial; weights = unit counts.
Estimator = maximum-weight compatible set (MWIS on the conflict graph, solved exactly with CP-SAT).

  part A  quartet phase diagram: recovery vs (long pendant rate p, internal rate r), exact boundary b = a
  part B  n = 10, three regimes: global separation (Thm 1), local separation (Thm 2: every false split
          lighter than every true split it crosses; must imply exact recovery in every replicate), and the
          pairwise Chernoff bound P(fail) <= sum_{tau x sigma} (1 - (sqrt m_sigma - sqrt m_tau)^2)^L (Thm 3),
          with m(.) estimated on 2e6 units; population optimum = MWIS on m itself (gamma = -1 flags a true tree that is not optimal); and the margin gamma
          (largest g keeping the true tree optimal after subtracting g from every true split) with the
          uniform Bernstein bound P(fail) <= sum_sigma 2 exp(-L gamma^2 / (8 m_sigma + 4 gamma / 3)) (Thm 4)
Usage: sim_identifiability.py OUT_DIR
"""
import sys, os, json, itertools, math
import numpy as np
from ortools.sat.python import cp_model

OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261118)


def canon(Y, n):
    """unrooted split of {0..n-1} with 1-set Y, as the frozenset side not containing witness 0"""
    Y = frozenset(Y)
    return Y if 0 not in Y else frozenset(range(n)) - Y


def nontrivial(Y, n): return 2 <= len(Y) <= n - 2


def compatible(A, B, n):
    U = frozenset(range(n)); Ac, Bc = U - A, U - B
    return not (A & B and A & Bc and Ac & B and Ac & Bc)


def mwis(weights, n):
    S = list(weights)
    if not S: return frozenset()
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in S]
    for i, j in itertools.combinations(range(len(S)), 2):
        if not compatible(S[i], S[j], n): m.Add(x[i] + x[j] <= 1)
    m.Maximize(sum(weights[s] * x[i] for i, s in enumerate(S)))
    sol = cp_model.CpSolver(); sol.parameters.num_workers = 8; st = sol.Solve(m)
    assert st == cp_model.OPTIMAL
    return frozenset(s for i, s in enumerate(S) if sol.Value(x[i]))


# ---------- rooted tree: list of edges (child clade as frozenset of leaves) ----------
class Tree:
    def __init__(self, clades, n, rates):
        self.n, self.clades, self.rates = n, clades, np.array(rates)
        self.M = np.array([[1 if l in c else 0 for l in range(n)] for c in clades], dtype=bool)  # edges x leaves
        self.true = {canon(c, n) for c in clades if nontrivial(c, n)}

    def units(self, L, eta=0.0):
        inn = rng.random((L, len(self.clades))) < self.rates            # innovations per unit per edge
        Y = (inn[:, :, None] & self.M[None]).any(1)                      # leaf has reading 1
        out = [frozenset(np.flatnonzero(y)) for y in Y]
        if eta:
            tc = [c for c in self.clades if nontrivial(c, self.n)]
            for u in np.flatnonzero(rng.random(L) < eta):
                c = set(tc[rng.integers(len(tc))]); w = int(rng.integers(self.n))
                c ^= {w}; out[u] = frozenset(c)
        return out


def weigh(units, n):
    w = {}
    for Y in units:
        if nontrivial(Y, n):
            s = canon(Y, n); w[s] = w.get(s, 0) + 1
    return w


# ---------- part A: quartet ----------
def quartet_exact(p, ps, r):
    """exact per-unit probabilities of the three quartet splits; root on the internal edge (halves r, r)"""
    # edges: pendant 0,1,2,3 (rates p, ps, p, ps: long branches on 0 and 2), r1 above {0,1}, r2 above {2,3}
    rates = [p, ps, p, ps, r, r]; clades = [{0}, {1}, {2}, {3}, {0, 1}, {2, 3}]
    m = {}
    for pat in itertools.product([0, 1], repeat=6):
        pr = math.prod(q if b else 1 - q for q, b in zip(rates, pat))
        Y = set().union(*[c for c, b in zip(clades, pat) if b]) if any(pat) else set()
        if len(Y) == 2:
            s = canon(Y, 4); m[s] = m.get(s, 0) + pr
    return m


def part_a():
    ps, L, R = 0.01, 200, 300
    T = frozenset({1}); F = frozenset({2}); G = frozenset({3})   # canon sides not containing 0: 01|23 -> {2,3}; check below
    res = []
    for p in np.linspace(0.02, 0.6, 15):
        for r in np.geomspace(0.002, 0.1, 15):
            m = quartet_exact(p, ps, r)
            true_s = canon({0, 1}, 4); false_s = canon({0, 2}, 4)
            b, a = m.get(true_s, 0), max(m.get(false_s, 0), m.get(canon({0, 3}, 4), 0))
            tree = Tree([frozenset(c) for c in [{0}, {1}, {2}, {3}, {0, 1}, {2, 3}]], 4, [p, ps, p, ps, r, r])
            ok = 0
            for _ in range(R):
                w = weigh(tree.units(L), 4)
                best = max(w.items(), key=lambda kv: (kv[1], rng.random()))[0] if w else None
                ok += best == true_s
            res.append(dict(p=float(p), r=float(r), b=b, a=a, recovery=ok / R))
    json.dump(dict(ps=ps, L=L, reps=R, grid=res), open(os.path.join(OUT, 'partA_quartet.json'), 'w'), indent=1)
    # does the exact boundary b = a separate recovery > 1/2 from recovery < 1/2?
    agree = sum((g['b'] > g['a']) == (g['recovery'] > 0.5) for g in res)
    print(f"part A: boundary b>a agrees with recovery>1/2 at {agree}/{len(res)} grid points")


# ---------- part B: n = 10 ----------
def balanced10():
    # ((((0,1),2),(3,4)),(((5,6),7),(8,(9)))) rooted; clades of every non-root vertex
    cl = [{i} for i in range(10)] + [{0, 1}, {0, 1, 2}, {3, 4}, {0, 1, 2, 3, 4}, {5, 6}, {5, 6, 7}, {8, 9}, {5, 6, 7, 8, 9}]
    return [frozenset(c) for c in cl]


def separation(w, true, n):
    """(global, local): every true > every false; every false < every true split it crosses"""
    tw = {s: w.get(s, 0) for s in true}; fw = {s: v for s, v in w.items() if s not in true}
    glob = min(tw.values()) > max(fw.values(), default=0)
    loc = all(v < tw[t] for s, v in fw.items() for t in true if not compatible(s, t, n))
    return glob, loc


def gamma_margin(m, true, n, N):
    """largest g such that the true split set stays a maximum after subtracting g from each true split"""
    val = lambda S, g: sum(m.get(x, 0) - (g if x in true else 0) for x in S)
    def ok(g):
        w = {x: round((v - (g if x in true else 0)) * N) for x, v in m.items()}
        w = {x: v for x, v in w.items() if v > 0}
        S = mwis(w, n)
        return val(S, g) <= val(true, g) + 1 / N
    lo, hi = 0.0, min(m.get(x, 0) for x in true)
    if not ok(1e-9): return -1.0
    for _ in range(20):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def regime(name, rates, eta, Ls, R=100):
    n = 10; clades = balanced10(); tree = Tree(clades, n, rates)
    N = 2_000_000; m = {s: c / N for s, c in weigh(tree.units(N, eta), n).items()}
    b = min(m.get(s, 0) for s in tree.true); a = max(v for s, v in m.items() if s not in tree.true)
    pairs = [(m[s], m.get(t, 0)) for s in m if s not in tree.true for t in tree.true if not compatible(s, t, n)]
    lmargin = min(mt - ms for ms, mt in pairs)
    pop_opt = mwis({s: round(v * N) for s, v in m.items()}, n) == tree.true
    gam = gamma_margin(m, tree.true, n, N)
    print(f"{name}: b={b:.4f} a={a:.4f} local margin={lmargin:+.4f} population optimum = true tree: {pop_opt} gamma={gam:.4f}")
    rows = []
    for L in Ls:
        rec = gs = ls = ls_wrong = 0
        for _ in range(R):
            w = weigh(tree.units(L, eta), n); g, l = separation(w, tree.true, n)
            got = mwis(w, n) == tree.true
            rec += got; gs += g; ls += l; ls_wrong += l and not got
        bound = min(1.0, sum((1 - (math.sqrt(mt) - math.sqrt(ms)) ** 2) ** L for ms, mt in pairs if mt > ms)
                    + sum(1 for ms, mt in pairs if mt <= ms)
                    + sum((1 - m.get(t, 0)) ** L for t in tree.true))          # a true split with no support
        bern = min(1.0, sum(2 * math.exp(-L * gam ** 2 / (8 * v + 4 * gam / 3)) for v in m.values())) if gam > 0 else 1.0
        rows.append(dict(L=L, reps=R, recovery=rec / R, global_sep=gs / R, local_sep=ls / R,
                         local_sep_but_wrong=ls_wrong, fail_bound_pairwise=bound, fail_bound_gamma=bern))
        print(f"  L={L:5d} recovery={rec/R:.3f} global={gs/R:.3f} local={ls/R:.3f} local&wrong={ls_wrong} "
              f"bound pairwise {bound:.3g} gamma {bern:.3g} obs fail {1-rec/R:.3f}")
    return dict(name=name, rates=rates, eta=eta, b=b, a=a, local_margin=lmargin, gamma=gam,
                true_splits=['|'.join(map(str, sorted(x))) for x in tree.true],
                m={'|'.join(map(str, sorted(x))): v for x, v in m.items()}, population_optimum_true=pop_opt,
                n_conflicting_pairs=len(pairs), rows=rows)


def part_b():
    clades = balanced10(); pend = lambda c: len(c) == 1
    base = [0.10 if pend(c) else 0.04 for c in clades]
    lba = list(base); lba[0] = lba[3] = 0.45                        # two heavily corrupted witnesses, far apart
    hard = [0.15 if pend(c) else 0.03 for c in clades]; hard[0] = hard[3] = 0.3
    out = [regime('clean', base, 0.02, [100, 400, 1600, 6400]),
           regime('noisy', hard, 0.05, [800, 3200, 12800]),
           regime('long-branch', lba, 0.02, [800, 3200])]
    json.dump(dict(clades=[sorted(c) for c in clades], regimes=out), open(os.path.join(OUT, 'partB_n10.json'), 'w'), indent=1)


if __name__ == '__main__':
    part_a(); part_b()
