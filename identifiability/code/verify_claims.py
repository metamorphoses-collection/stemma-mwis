"""Machine checks of the combinatorial claims of the identifiability paper. Every check can fail.
Writes data/verify_claims.json and prints one line per claim.

  C1  Universality: the construction realises random graphs H exactly as conflict graphs (distinct splits).
  C2  Pairs exact, triples not: exhaustive on n = 4, 5, 6 witnesses (pairs completable; no obstructing triple
      on 4; count of obstructing triples on 5).
  C3  Displacement matching: random binary trees and random compatible sets S; |B| <= |A|, Hall matching
      exists, A nonempty for S != true.
  C4  Maximality: every non-trivial split compatible with all splits of a random binary tree is one of them;
      every compatible set has at most n-3 splits (random greedy maximal sets).
  C5  Rooted analogue: laminar families of subsets of size >= 2 have at most n-1 members, and the cluster
      system of a binary rooted tree is maximal.
  C6  Quartet formulas of the irreversible model equal exhaustive enumeration (random rates).
  C7  Irreversible model, general n: m(true split of e) >= pi_e prod_{f outside subtree(e), f != e}(1-pi_f) and
      total false mass <= sum over incomparable edge pairs pi_e pi_f (exhaustive enumeration, n = 6, 7).
  C8  gamma: the bisection characterisation agrees with the ratio definition (brute force over all compatible
      sets, n = 6).
  C9  Cutting planes: on random weighted families of partial splits (n = 5, 6) the lazy-cut loop (master solved
      by brute force) returns the weight of a brute-force maximum jointly completable family.
  C10 Certificates as optimisation statements: random binary trees (n = 5..7), random weights satisfying the
      local (resp. global) hypothesis; the true split set is the unique maximum-weight compatible set
      (brute force over all compatible sets).
  C11 Proper-partial obstructions: on five witnesses, obstructing triples of proper partial splits (quartets).
  C12 The full conflict graph on five witnesses is isomorphic to the line graph of K5 (complement of the Petersen
      graph); it has an induced 4-cycle and an induced 5-cycle (not chordal, not perfect). Full conflict graphs
      on n = 4..7 witnesses: vertex counts 2^(n-1)-n-1.
  C13 Crossing set: for random binary trees and every false split tau = A|B, the true splits crossed by tau are
      the internal edges of T[A] cap T[B], non-empty and connected.
  C14 Dynamic programme over laminar families (O(3^n)) returns the brute-force maximum-weight compatible set
      value on random weights, n = 5..8.
These are finite checks on random or exhaustive instances, not proofs. The script exits non-zero on any failure.
Usage: verify_claims.py OUT_DIR
"""
import sys, os, json, itertools, random, math
import networkx as nx
from ortools.sat.python import cp_model

OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
random.seed(20261118)
R = {}


def compatible(A, B, X):
    Ac, Bc = X - A, X - B
    return not (A & B and A & Bc and Ac & B and Ac & Bc)


def canon(A, X):
    A = frozenset(A); return A if min(X) not in A else frozenset(X) - A


def nontrivial(A, X): return 2 <= len(A) <= len(X) - 2


# ---------- C1 ----------
def universality(H):
    V, E = list(H.nodes), [tuple(e) for e in H.edges]
    W = [('p', v, i) for v in V for i in (0, 1)] + [('q', e, a, b) for e in E for a in (0, 1) for b in (0, 1)] + [('z', 1), ('z', 2)]
    X = frozenset(W)
    sides = {}
    for v in V:
        s = {('p', v, 0), ('p', v, 1)}
        for e in E:
            if v in e:
                k = e.index(v)
                s |= {('q', e, a, b) for a in (0, 1) for b in (0, 1) if (a, b)[k] == 1}
        sides[v] = frozenset(s)
    ok = len(X) == 4 * len(E) + 2 * len(V) + 2
    ok &= len({canon(s, X) for s in sides.values()}) == len(V)
    ok &= all(nontrivial(s, X) for s in sides.values())
    for u, v in itertools.combinations(V, 2):
        ok &= (not compatible(sides[u], sides[v], X)) == H.has_edge(u, v)
    return ok


c1 = []
for _ in range(300):
    nv = random.randint(2, 8); H = nx.gnp_random_graph(nv, random.random(), seed=random.randint(0, 10**9))
    c1.append(universality(H))
R['C1_universality'] = dict(graphs=len(c1), all_realised=all(c1))


# ---------- C2 ----------
def partial_splits(n):
    out = []
    X = list(range(n))
    for k in range(4, n + 1):
        for Y in itertools.combinations(X, k):
            Y = frozenset(Y)
            for r in range(2, k - 1):
                for A in itertools.combinations(sorted(Y), r):
                    A = frozenset(A)
                    if min(Y) in A: out.append((A, Y - A, Y))
    return out


def overlap_compat(s, t):
    Z = s[2] & t[2]
    return compatible(s[0] & Z, t[0] & Z, Z)


def completable(fam, n):
    X = frozenset(range(n))
    choices = []
    for A, B, Y in fam:
        un = sorted(X - Y)
        choices.append([A | frozenset(c) for r in range(len(un) + 1) for c in itertools.combinations(un, r)])
    for ext in itertools.product(*choices):
        if all(compatible(a, b, X) for a, b in itertools.combinations(ext, 2)): return True
    return False


c2 = {}
for n in (4, 5, 6):
    P = partial_splits(n)
    pairs_ok = all(completable([s, t], n) for s, t in itertools.combinations(P, 2) if overlap_compat(s, t))
    obstr = None
    if n <= 5:
        obstr = sum(1 for T in itertools.combinations(P, 3)
                    if all(overlap_compat(a, b) for a, b in itertools.combinations(T, 2)) and not completable(list(T), n))
    c2[n] = dict(partial_splits=len(P), pairs_completable=pairs_ok, obstructing_triples=obstr)
ex = [(frozenset({0, 1}), frozenset({2, 3}), frozenset({0, 1, 2, 3})), (frozenset({0, 1}), frozenset({2, 4}), frozenset({0, 1, 2, 4})),
      (frozenset({0, 3}), frozenset({1, 4}), frozenset({0, 1, 3, 4}))]
c2['example_pqr'] = dict(pairwise=all(overlap_compat(a, b) for a, b in itertools.combinations(ex, 2)), completable=completable(ex, 5))
R['C2_pairs_triples'] = c2


# ---------- trees ----------
def random_binary_tree(n):
    """unrooted binary tree as list of non-trivial splits (canonical) built by random leaf insertion"""
    G = nx.Graph(); G.add_edges_from([(0, 'h'), (1, 'h'), (2, 'h')]); k = 0
    for leaf in range(3, n):
        u, v = random.choice(list(G.edges)); w = f'i{k}'; k += 1
        G.remove_edge(u, v); G.add_edges_from([(u, w), (w, v), (w, leaf)])
    X = frozenset(range(n)); S = set()
    for u, v in G.edges:
        H = G.copy(); H.remove_edge(u, v)
        side = frozenset(x for x in nx.node_connected_component(H, u) if isinstance(x, int))
        if nontrivial(side, X): S.add(canon(side, X))
    return S


def all_splits(n):
    X = frozenset(range(n))
    return [frozenset(A) for r in range(2, n - 1) for A in itertools.combinations(range(1, n), r) if nontrivial(frozenset(A), X)]


def random_compatible(cands, X, p):
    S = []
    for s in random.sample(cands, len(cands)):
        if random.random() < p and all(compatible(s, t, X) for t in S): S.append(s)
    return S


c3 = dict(trials=0, size_ok=True, A_nonempty=True, matching=True)
c4 = dict(trees=0, maximal=True, max_size_ok=True)
for _ in range(400):
    n = random.randint(4, 9); X = frozenset(range(n)); T = random_binary_tree(n); C = all_splits(n)
    c4['trees'] += 1
    c4['maximal'] &= len(T) == n - 3 and all(s in T for s in C if all(compatible(s, t, X) for t in T))
    for _ in range(10):
        S = set(random_compatible(C, X, random.random()))
        c4['max_size_ok'] &= len(S) <= n - 3
        if S == T: continue
        A, B = T - S, S - T; c3['trials'] += 1
        c3['size_ok'] &= len(B) <= len(A); c3['A_nonempty'] &= bool(A)
        G = nx.Graph(); G.add_nodes_from([('b', x) for x in B]); G.add_nodes_from([('a', x) for x in A])
        G.add_edges_from((('b', x), ('a', y)) for x in B for y in A if not compatible(x, y, X))
        M = nx.bipartite.maximum_matching(G, top_nodes=[('b', x) for x in B]) if B else {}
        c3['matching'] &= sum(1 for k in M if k[0] == 'b') == len(B)
R['C3_displacement_matching'] = c3; R['C4_maximality'] = c4


# ---------- C5 rooted ----------
def laminar(A, B): return not (A & B) or A <= B or B <= A


def random_rooted_clusters(n):
    nodes = [frozenset({i}) for i in range(n)]; cl = set()
    while len(nodes) > 1:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; cl.add(c); nodes.append(c)
    return cl


c5 = dict(trees=0, size_ok=True, maximal=True, laminar_bound=True)
for _ in range(300):
    n = random.randint(3, 8); Cl = random_rooted_clusters(n)
    subsets = [frozenset(A) for r in range(2, n + 1) for A in itertools.combinations(range(n), r)]
    c5['trees'] += 1; c5['size_ok'] &= len(Cl) == n - 1
    c5['maximal'] &= all(s in Cl for s in subsets if all(laminar(s, c) for c in Cl))
    F = []
    for s in random.sample(subsets, len(subsets)):
        if all(laminar(s, t) for t in F): F.append(s)
    c5['laminar_bound'] &= len(F) <= n - 1
R['C5_rooted'] = c5


# ---------- irreversible model, exact enumeration ----------
def enumerate_m(clades, rates, n):
    X = frozenset(range(n)); m = {}
    for pat in itertools.product([0, 1], repeat=len(clades)):
        pr = math.prod(q if b else 1 - q for q, b in zip(rates, pat))
        Y = frozenset().union(*[c for c, b in zip(clades, pat) if b]) if any(pat) else frozenset()
        if nontrivial(Y, X): s = canon(Y, X); m[s] = m.get(s, 0) + pr
    return m


c6 = []
for _ in range(200):
    p = [random.random() * 0.6 for _ in range(4)]; r1, r2 = random.random() * 0.3, random.random() * 0.3
    clades = [frozenset(c) for c in [{0}, {1}, {2}, {3}, {0, 1}, {2, 3}]]
    m = enumerate_m(clades, p + [r1, r2], 4)
    b = lambda x: 1 - x
    t = b(r2) * b(p[2]) * b(p[3]) * (r1 + b(r1) * p[0] * p[1]) + b(r1) * b(p[0]) * b(p[1]) * (r2 + b(r2) * p[2] * p[3])
    f13 = b(r1) * b(r2) * (p[0] * p[2] * b(p[1]) * b(p[3]) + p[1] * p[3] * b(p[0]) * b(p[2]))
    f14 = b(r1) * b(r2) * (p[0] * p[3] * b(p[1]) * b(p[2]) + p[1] * p[2] * b(p[0]) * b(p[3]))
    X = frozenset(range(4))
    c6.append(abs(m.get(canon({0, 1}, X), 0) - t) < 1e-12 and abs(m.get(canon({0, 2}, X), 0) - f13) < 1e-12
              and abs(m.get(canon({0, 3}, X), 0) - f14) < 1e-12)
R['C6_quartet_formulas'] = dict(trials=len(c6), all_equal=all(c6))


def random_rooted_edges(n):
    """clades below every edge of a random rooted binary tree (root edge excluded: root is a vertex)"""
    cl = random_rooted_clusters(n)
    return [frozenset({i}) for i in range(n)] + [c for c in cl if len(c) < n]


c7 = dict(trials=0, true_lower=True, false_upper=True)
for _ in range(60):
    n = random.choice([6, 7]); clades = random_rooted_edges(n); X = frozenset(range(n))
    rates = [random.random() * 0.4 for _ in clades]
    m = enumerate_m(clades, rates, n)
    true = {canon(c, X) for c in clades if nontrivial(c, X)}
    c7['trials'] += 1
    for c, pe in zip(clades, rates):
        if not nontrivial(c, X): continue
        lower = pe * math.prod(1 - q for d, q in zip(clades, rates) if not d <= c)
        c7['true_lower'] &= m.get(canon(c, X), 0) >= lower - 1e-12
    incomparable = sum(qe * qf for (e, qe), (f, qf) in itertools.combinations(zip(clades, rates), 2) if not (e & f))
    c7['false_upper'] &= sum(v for s, v in m.items() if s not in true) <= incomparable + 1e-12
R['C7_irreversible_bounds'] = c7


# ---------- C8 gamma ----------
def compatible_sets(C, X):
    out = [frozenset()]
    def rec(i, cur):
        for j in range(i, len(C)):
            if all(compatible(C[j], t, X) for t in cur):
                nxt = cur + [C[j]]; out.append(frozenset(nxt)); rec(j + 1, nxt)
    rec(0, []); return out


c8 = dict(trials=0, agree=True, max_abs_diff=0.0, tolerance=1e-9)
for _ in range(25):
    n = 6; X = frozenset(range(n)); T = random_binary_tree(n); C = all_splits(n)
    m = {s: random.random() * (1.0 if s in T else 0.6) for s in C}
    Ss = [S for S in compatible_sets(C, X) if S != T]
    ratio = min((sum(m[s] for s in T - S) - sum(m[s] for s in S - T)) / len(T - S) for S in Ss)
    def optimal(g): return all(sum(m[s] - (g if s in T else 0) for s in T) >= sum(m[s] - (g if s in T else 0) for s in S) - 1e-12 for S in Ss)
    lo, hi = -2.0, 2.0
    for _ in range(50):
        mid = (lo + hi) / 2; lo, hi = (mid, hi) if optimal(mid) else (lo, mid)
    c8['trials'] += 1; d = abs(lo - ratio); c8['max_abs_diff'] = max(c8['max_abs_diff'], d); c8['agree'] &= d < 1e-9
R['C8_gamma'] = c8


# ---------- C9 cutting planes ----------
def brute_best(fam, w, n, ok):
    best, bestS = -1, None
    for r in range(len(fam) + 1):
        for S in itertools.combinations(range(len(fam)), r):
            if ok(S):
                v = sum(w[i] for i in S)
                if v > best: best, bestS = v, S
    return best, bestS


c9 = dict(instances=0, agree=True, cuts_total=0)
for _ in range(60):
    n = random.choice([5, 6]); P = partial_splits(n); fam = random.sample(P, min(len(P), random.randint(5, 9)))
    w = [random.randint(1, 9) for _ in fam]
    exact, _ = brute_best(fam, w, n, lambda S: completable([fam[i] for i in S], n))
    cuts = []
    while True:
        pair_ok = lambda S: all(overlap_compat(fam[i], fam[j]) for i, j in itertools.combinations(S, 2)) and \
            all(not set(K) <= set(S) for K in cuts)
        val, S = brute_best(fam, w, n, pair_ok)
        if completable([fam[i] for i in S], n): break
        cuts.append(S)
    c9['instances'] += 1; c9['cuts_total'] += len(cuts); c9['agree'] &= val == exact
R['C9_cutting_planes'] = c9


# ---------- C10 certificates ----------
c10 = dict(instances_local=0, instances_global=0, unique_true_optimum=True)
for _ in range(80):
    n = random.choice([5, 6, 7]); X = frozenset(range(n)); T = random_binary_tree(n); C = all_splits(n)
    for kind in ('local', 'global'):
        w = {}
        base = {s: random.randint(5, 12) for s in T}
        for s in C:
            if s in T: w[s] = base[s]
            else:
                cap = min(base[t] for t in T if not compatible(s, t, X)) if kind == 'local' else min(base.values())
                w[s] = random.randint(0, cap - 1)
        Ss = compatible_sets(C, X)
        vals = sorted(((sum(w[x] for x in S), S) for S in Ss), key=lambda t: -t[0])
        c10['instances_' + kind] += 1
        c10['unique_true_optimum'] &= vals[0][1] == frozenset(T) and vals[1][0] < vals[0][0]
R['C10_certificates'] = c10


# ---------- C11 proper-partial obstructions on five witnesses ----------
P5 = [s for s in partial_splits(5) if len(s[2]) < 5]
R['C11_proper_partial_triples_n5'] = sum(1 for T in itertools.combinations(P5, 3)
    if all(overlap_compat(a, b) for a, b in itertools.combinations(T, 2)) and not completable(list(T), 5))


# ---------- C12 full conflict graphs ----------
def full_conflict_graph(n):
    X = frozenset(range(n)); C = all_splits(n); G = nx.Graph(); G.add_nodes_from(C)
    G.add_edges_from((a, b) for a, b in itertools.combinations(C, 2) if not compatible(a, b, X)); return G


def has_induced_cycle(G, k):
    for cyc in itertools.permutations(G.nodes, k):
        if cyc[0] != min(cyc, key=lambda s: sorted(s)): continue
        H = G.subgraph(cyc)
        if H.number_of_edges() == k and all(H.degree(v) == 2 for v in cyc) and nx.is_connected(H): return True
    return False


G5 = full_conflict_graph(5)
R['C12_structure'] = dict(n5_iso_line_graph_K5=nx.is_isomorphic(G5, nx.line_graph(nx.complete_graph(5))),
                          n5_iso_complement_petersen=nx.is_isomorphic(G5, nx.complement(nx.petersen_graph())),
                          n5_chordal=nx.is_chordal(G5), n5_induced_C4=has_induced_cycle(G5, 4), n5_induced_C5=has_induced_cycle(G5, 5),
                          vertex_counts={n: full_conflict_graph(n).number_of_nodes() for n in range(4, 8)})


# ---------- C13 crossing set ----------
def tree_graph(n):
    G = nx.Graph(); G.add_edges_from([(0, 'h'), (1, 'h'), (2, 'h')]); k = 0
    for leaf in range(3, n):
        u, v = random.choice(list(G.edges)); w_ = f'i{k}'; k += 1
        G.remove_edge(u, v); G.add_edges_from([(u, w_), (w_, v), (w_, leaf)])
    return G


def edge_split(G, e, X):
    H = G.copy(); H.remove_edge(*e)
    return frozenset(x for x in nx.node_connected_component(H, e[0]) if isinstance(x, int))


def span_edges(G, A):
    A = list(A); ed = set()
    for a, b in itertools.combinations(A, 2):
        path = nx.shortest_path(G, a, b); ed |= {frozenset(e) for e in zip(path, path[1:])}
    return ed


c13 = dict(trees=0, false_splits=0, crossing_equals_intersection=True, nonempty=True, connected=True)
for _ in range(120):
    n = random.randint(5, 9); X = frozenset(range(n)); G = tree_graph(n)
    internal = [e for e in G.edges if not (isinstance(e[0], int) or isinstance(e[1], int))]
    T = {canon(edge_split(G, e, X), X): frozenset(e) for e in internal}
    c13['trees'] += 1
    for t in all_splits(n):
        if t in T: continue
        A, B = t, X - t
        crossed = {T[s] for s in T if not compatible(t, s, X)}
        inter = span_edges(G, A) & span_edges(G, B)
        c13['false_splits'] += 1
        c13['crossing_equals_intersection'] &= crossed == inter
        c13['nonempty'] &= bool(crossed)
        H = nx.Graph(); H.add_edges_from(tuple(e) for e in crossed)
        c13['connected'] &= nx.is_connected(H)
R['C13_crossing_set'] = c13


# ---------- C14 dynamic programme ----------
def dp_best(n, w):
    # root at witness 0; compatible non-trivial splits <-> laminar families of sides avoiding 0, sizes 2..n-2
    Y = list(range(1, n)); idx = {y: i for i, y in enumerate(Y)}; full = (1 << len(Y)) - 1
    wt = {}
    for s, v in w.items():
        side = s if 0 not in s else frozenset(range(n)) - s
        wt[sum(1 << idx[y] for y in side)] = v
    F, P = {0: 0}, {0: 0}
    for C in sorted(range(1, full + 1), key=lambda c: bin(c).count('1')):
        low = C & -C; rest = C ^ low
        best_inner = 0; sub = rest
        while True:                      # D = low | sub, proper subset of C
            D = low | sub
            if D != C: best_inner = max(best_inner, F[D] + P[C ^ D])
            if sub == 0: break
            sub = (sub - 1) & rest
        F[C] = best_inner + max(0, wt.get(C, 0))
        bp = 0; sub = rest
        while True:
            D = low | sub; bp = max(bp, F[D] + P[C ^ D])
            if sub == 0: break
            sub = (sub - 1) & rest
        P[C] = bp
    return F[full]


c14 = dict(instances=0, agree=True)
for _ in range(40):
    n = random.randint(5, 8); X = frozenset(range(n)); C = all_splits(n)
    w = {s: random.choice([0, 0, random.randint(1, 20)]) for s in C}
    G = nx.Graph(); G.add_nodes_from([s for s in C if w[s] > 0])
    G.add_edges_from((a, b) for a, b in itertools.combinations(list(G.nodes), 2) if not compatible(a, b, X))
    m = cp_model.CpModel(); x = {v: m.NewBoolVar('') for v in G.nodes}
    for a, b in G.edges: m.Add(x[a] + x[b] <= 1)
    m.Maximize(sum(w[v] * x[v] for v in G.nodes)); sol = cp_model.CpSolver(); sol.Solve(m)
    c14['instances'] += 1; c14['agree'] &= dp_best(n, w) == round(sol.ObjectiveValue())
R['C14_dynamic_programme'] = c14

json.dump(R, open(os.path.join(OUT, 'verify_claims.json'), 'w'), indent=1, default=str)
for k, v in R.items(): print(k, v)
expect = [R['C1_universality']['all_realised'], R['C2_pairs_triples'][4]['obstructing_triples'] == 0,
          R['C2_pairs_triples'][5]['obstructing_triples'] == 200, all(R['C2_pairs_triples'][k]['pairs_completable'] for k in (4, 5, 6)),
          R['C2_pairs_triples']['example_pqr'] == dict(pairwise=True, completable=False),
          all(v for k, v in R['C3_displacement_matching'].items() if k != 'trials'),
          R['C4_maximality']['maximal'] and R['C4_maximality']['max_size_ok'],
          all(v for k, v in R['C5_rooted'].items() if k != 'trees'), R['C6_quartet_formulas']['all_equal'],
          R['C7_irreversible_bounds']['true_lower'] and R['C7_irreversible_bounds']['false_upper'],
          R['C8_gamma']['agree'], R['C9_cutting_planes']['agree'], R['C10_certificates']['unique_true_optimum'],
          R['C11_proper_partial_triples_n5'] == 140,
          R['C12_structure']['n5_iso_complement_petersen'] and not R['C12_structure']['n5_chordal']
          and R['C12_structure']['n5_induced_C4'] and R['C12_structure']['n5_induced_C5']
          and R['C12_structure']['vertex_counts'] == {n: 2 ** (n - 1) - n - 1 for n in range(4, 8)},
          all(v for k, v in R['C13_crossing_set'].items() if k not in ('trees', 'false_splits')),
          R['C14_dynamic_programme']['agree']]
print('ALL PASS' if all(expect) else f'FAIL at {[i for i, e in enumerate(expect) if not e]}')
sys.exit(0 if all(expect) else 1)
