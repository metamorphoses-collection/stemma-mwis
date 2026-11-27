r"""Does the Bedier effect reproduce the excess resolution of editors' stemmata?

Traditions are simulated with the birth-death null model of Camps, Randon-Furling & Godreau (PNAS Nexus 5(7) pgag207,
2026; preprint arXiv 2505.19246): one original at t = 0, each living manuscript copied at rate lam and lost at rate mu
during an active phase of T_a, then lost at rate mu during a decimation phase of T_d (no copying), simulated exactly in
continuous time (Gillespie; a wait that crosses T_a restarts at T_a with the decimation rate mu*n, by memorylessness).
PARAMS=published (default): the published constant-rate SBI estimates (Table 1), lam = 0.032, mu = 0.014 per year,
T_a = T_d = 250 years. PARAMS=preprint: the earlier setting, lam = 8.5e-3, mu = 3.8e-3 per three-month step,
T_a = T_d = 1000 steps (midpoints of the preprint's ranges; kept for reproducing the old outputs).
Traditions with 3 to 30 surviving manuscripts are kept. Every CP-SAT solve must return OPTIMAL; a tradition with any
non-optimal solve is excluded from all arms and counted in the output (non_optimal_excluded).

The true stemma is the full tree pruned to the survivors: lost leaves removed recursively and lost nodes with one
child contracted. A surviving manuscript with descendants is represented as a leaf hanging from its node.

Texts: L word positions. Every copy changes each position with probability e; a change takes the position's
attractive reading with probability alpha and otherwise a reading of its own (attraction model). The archetype
text is the reference (rooted coding): a unit proposes the cluster of surviving manuscripts that share a departure
from it, of size 2..n-1. Invariant positions propose nothing.

Reconstructions: (raw) maximum-weight laminar family of clusters, exact (CP-SAT); (rival) the same, then every
selected cluster must beat its strongest crossing candidate by the exact one-sided binomial test at level 0.05
(Proposition rival), a crossing weight 0 if there is none; clusters failing it are removed (which keeps the family
laminar); (k>=K) the Bedier threshold of Theorem bedier with c = K/L: the same optimum restricted to clusters with
at least K supporting units, K = 2, 3, 5.

Metrics on every rooted tree, true or reconstructed, in the same leaf representation: bifid root (root with two
children), shares of internal nodes with 2, exactly 3 and >= 3 children, and i3 = share of leaf triples that some cluster resolves
(contains exactly two of them), the index of Camps et al. Usage: bedier_on_birthdeath.py OUT.json [n_traditions]. Environment: E (mean per-copy change rate, default 0.02),
L (positions, default 1000), SHAPE (gamma shape of per-position rates; 0 = equal rates, default; Spencer et al. 2004
measured 0.38 on an artificial tradition), ALPHAS (comma list, default 0.3,0.6,0.9)."""
import os
import sys, json, math, random, itertools
import numpy as np
from ortools.sat.python import cp_model

OUT = sys.argv[1]; NT = int(sys.argv[2]) if len(sys.argv) > 2 else 300
PARAMS = os.environ.get('PARAMS', 'published')
LAM, MU, TA, TD = {'published': (0.032, 0.014, 250, 250), 'preprint': (8.5e-3, 3.8e-3, 1000, 1000)}[PARAMS]
L, E = int(os.environ.get('L', 1000)), float(os.environ.get('E', 0.02))
SHAPE = float(os.environ.get('SHAPE', 0))
ALPHAS = [float(a) for a in os.environ.get('ALPHAS', '0.3,0.6,0.9').split(',')]
rng = random.Random(20261118); nrng = np.random.default_rng(20261118)


def simulate_tree():
    """Gillespie birth-death; returns parent list, alive flags (at the end)."""
    parent, alive = [-1], [True]
    living = [0]; t = 0.0
    while living:
        n = len(living)
        active = t < TA
        rate = (LAM + MU) * n if active else MU * n
        t += rng.expovariate(rate)
        if active and t >= TA: t = TA; continue                 # wait crossed the phase boundary: restart at TA
        if t >= TA + TD: break
        if active and rng.random() < LAM / (LAM + MU):
            p = rng.choice(living); parent.append(p); alive.append(True); living.append(len(parent) - 1)
        else:
            k = rng.randrange(n); d = living[k]; living[k] = living[-1]; living.pop(); alive[d] = False
        if len(parent) > 5000: return None
    return parent, alive


def prune(parent, alive):
    """Clusters (sets of survivor indices) of the pruned rooted stemma, with survivors as leaves."""
    n = len(parent); children = [[] for _ in range(n)]
    for i, p in enumerate(parent):
        if p >= 0: children[p].append(i)
    surv = [i for i in range(n) if alive[i]]; sidx = {s: k for k, s in enumerate(surv)}
    below = [None] * n
    for i in reversed(range(n)):                               # children have larger ids than parents
        s = {sidx[i]} if alive[i] else set()
        for c in children[i]: s |= below[c]
        below[i] = frozenset(s)
    clusters = set()
    for i in range(n):
        b = below[i]
        if 2 <= len(b) <= len(surv) - 1:
            # node i defines a cluster if it is not a contracted unary lost node: count non-empty child groups
            groups = sum(1 for c in children[i] if below[c]) + (1 if alive[i] else 0)
            if groups >= 2: clusters.add(b)
    return surv, clusters


def texts(parent, alpha):
    n = len(parent); T = np.zeros((n, L), dtype=np.int64)
    rates = np.minimum(E * nrng.gamma(SHAPE, 1 / SHAPE, L), 1.0) if SHAPE > 0 else np.full(L, E)
    for i in range(1, n):
        T[i] = T[parent[i]]
        ch = nrng.random(L) < rates
        att = nrng.random(L) < alpha
        T[i][ch & att] = 1
        T[i][ch & ~att] = i + 1                                 # a reading of its own
    return T


def candidates(T, surv):
    X = len(surv); w = {}
    S = T[surv]
    for u in range(L):
        col = S[:, u]
        for v in set(col.tolist()):
            if v == 0: continue
            C = frozenset(np.flatnonzero(col == v).tolist())
            if 2 <= len(C) <= X - 1: w[C] = w.get(C, 0) + 1
    return w


def laminar(a, b): return not (a & b) or a <= b or b <= a


class NonOptimal(Exception): pass


def mwis(w):
    V = list(w)
    if not V: return set()
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in V]
    for i, j in itertools.combinations(range(len(V)), 2):
        if not laminar(V[i], V[j]): m.Add(x[i] + x[j] <= 1)
    m.Maximize(sum(w[v] * x[k] for k, v in enumerate(V)))
    s = cp_model.CpSolver(); s.parameters.num_workers = 8; s.parameters.max_time_in_seconds = 60
    st = s.Solve(m)
    if st != cp_model.OPTIMAL: raise NonOptimal(s.StatusName(st))
    return {v for k, v in enumerate(V) if s.Value(x[k])}


def binom_tail(k, n): return sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n


def rival_filter(S, w, level=0.05):
    keep = set()
    for C in S:
        r = max((w[t] for t in w if not laminar(C, t)), default=0)
        if binom_tail(w[C], w[C] + r) < level: keep.add(C)
    return keep


def metrics(clusters, n):
    X = frozenset(range(n)); cl = set(clusters) | {X}
    def children(C):
        subs = [D for D in cl if D < C]
        maxs = [D for D in subs if not any(D < E2 for E2 in subs)]
        covered = set().union(*maxs) if maxs else set()
        return len(maxs) + len(C - covered)
    degs = [children(C) for C in cl]
    trip = list(itertools.combinations(range(n), 3))
    res = sum(1 for t in trip if any(len(C & set(t)) == 2 for C in clusters)) / len(trip)
    return dict(bifid_root=children(X) == 2, d2=sum(d == 2 for d in degs), d3=sum(d == 3 for d in degs),
                ge3=sum(d >= 3 for d in degs), internal=len(degs), i3=res)


def published_degrees(parent, alive):
    """Out-degrees of the internal nodes of the stemma built as in Camps et al. (Methods, Fig. 7): dead leaves removed
    recursively, then non-branching dead nodes contracted; surviving witnesses stay where they are, so a surviving
    exemplar is an internal node (out-degree = its surviving lines, possibly 1). Root = topmost remaining node."""
    n = len(parent); children = [[] for _ in range(n)]
    for i, p in enumerate(parent):
        if p >= 0: children[p].append(i)
    keep = [False] * n
    for i in reversed(range(n)):
        keep[i] = alive[i] or any(keep[c] for c in children[i])
    def lines(i):                                              # kept children after contracting unary dead nodes
        out = []
        for c in children[i]:
            if not keep[c]: continue
            while not alive[c] and sum(keep[g] for g in children[c]) == 1:
                c = next(g for g in children[c] if keep[g])
            out.append(c)
        return out
    r = 0
    while not alive[r] and sum(keep[g] for g in children[r]) == 1: r = next(g for g in children[r] if keep[g])
    degs, stack = [], [r]
    while stack:
        i = stack.pop(); ch = lines(i)
        if ch: degs.append(len(ch)); stack.extend(ch)
    return dict(pub_bifid_root=len(lines(r)) == 2, pub_ge3=sum(d >= 3 for d in degs), pub_internal=len(degs))


def summarise(rows):
    tot = lambda k: sum(r[k] for r in rows)
    out = dict(traditions=len(rows), bifid_root=round(sum(r['bifid_root'] for r in rows) / len(rows), 3),
               p_deg2=round(tot('d2') / tot('internal'), 3), p_deg3=round(tot('d3') / tot('internal'), 3),
               p_deg_ge3=round(tot('ge3') / tot('internal'), 3), i3=round(float(np.mean([r['i3'] for r in rows])), 3))
    if 'pub_internal' in rows[0]:                              # truth only: the published stemma convention
        out.update(pub_bifid_root=round(sum(r['pub_bifid_root'] for r in rows) / len(rows), 3),
                   pub_p_deg_ge3=round(tot('pub_ge3') / tot('pub_internal'), 3))
    return out


results = {}
for alpha in ALPHAS:
    KS = (2, 3, 5); arms = ['raw', 'rival'] + [f'k>={k}' for k in KS]
    rows = {a: [] for a in ['true'] + arms}; acc = {a: [0, 0, 0] for a in arms}
    got = 0; tries = 0; nonopt = []
    while got < NT and tries < 200000:
        tries += 1
        sim = simulate_tree()
        if sim is None: continue
        parent, alive = sim
        ns = sum(alive)
        if not 3 <= ns <= 30: continue
        surv, truth = prune(parent, alive)
        T = texts(parent, alpha)
        w = candidates(T, surv)
        try:
            raw = mwis(w); rec = {'raw': raw, 'rival': rival_filter(raw, w)}
            for k in KS: rec[f'k>={k}'] = mwis({C: x for C, x in w.items() if x >= k})
        except NonOptimal as ex:
            nonopt.append(dict(try_=tries, survivors=ns, candidates=len(w), status=str(ex))); continue
        rows['true'].append(dict(metrics(truth, ns), **published_degrees(parent, alive)))
        for k, S in rec.items():
            rows[k].append(metrics(S, ns))
            acc[k][0] += len(S & truth); acc[k][1] += len(S - truth); acc[k][2] += len(truth)
        got += 1
    results[f'alpha={alpha}'] = dict({k: summarise(v) for k, v in rows.items()},
                                     recovery={k: dict(true_found=a[0], false_added=a[1], true_total=a[2]) for k, a in acc.items()},
                                     simulated=tries, non_optimal_excluded=len(nonopt), non_optimal=nonopt)
    print(f'alpha={alpha}', json.dumps(results[f'alpha={alpha}']), flush=True)
results['params'] = dict(params=PARAMS, clock='gillespie, phase boundary restarted at T_active', lam=LAM, mu=MU, T_active=TA, T_decimation=TD, L=L, e=E, gamma_shape=SHAPE, traditions_per_alpha=NT,
                         survivors='3..30', camps_table1=dict(real=dict(bifid=0.77, p_deg2=0.81, p_deg3=0.11, i3=0.93),
                                                               simulated=dict(bifid='0.60-0.64', p_deg2='0.71-0.78', p_deg3='0.14-0.15', i3='0.80-0.85')),
                         camps2026_table2=dict(empirical=dict(bifid=0.77, bifid_ci='0.65-0.86', i3=0.93, i3_ci='0.89-0.96',
                                                              p_deg_ge3=0.36, p_deg_ge3_ci='0.32-0.38'),
                                               constant_rate=dict(bifid=0.65, i3=0.83, p_deg_ge3=0.33)))
json.dump(results, open(OUT, 'w'), indent=1)
