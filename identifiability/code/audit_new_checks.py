r"""New machine checks for results of the technical report that had none (final audit, 22 Nov 2026).

Each check is a function that returns a dict with 'passed' (bool), 'instances' and 'result' (one line each, LaTeX-safe,
used in the audit table) and 'numbers' (the key quantities); every check can fail. Seeded, independent of the proofs:
supports come from exact enumeration of the unit law (attraction model, all innovation outcomes) or from a designed
law, and the estimators are evaluated by brute force over all compatible sets.

  exclusive       Lemma exclusive: splits (and clusters) proposed by one unit are pairwise compatible.
  nphard          Cor. nphard: the construction of Thm universal preserves MWIS values (and the optimal sets).
  buneman         Prop. buneman: split set -> tree -> split set round trip; uniqueness up to isomorphism.
  threshold       Prop. threshold: restricted optimum equals the true set (all compatible sets enumerated).
  pairwise        Thm pairwise: the bound >= simulated P(S_hat != T) (ties count as errors).
  consistency     Thm consistency: (a) both bounds >= simulated error; (b) P(S_hat = T) -> 0 when m(S) > m(T).
  samplesize      Cor. samplesize: the formula's L gives a Bernstein bound <= eps; simulated error <= eps.
  blocks          Prop. blocks: J i.i.d. blocks with dependence inside: the bound >= simulated error.
  gammablocks     Prop. gammablocks: simultaneous coverage of gamma_lo <= gamma >= 1 - delta; false certificates <= delta.
  textillegible   Prop. textillegible: reading-dependent illegibility; both bounds of Thm fragmentary(b) with delta,
                  exact enumeration; the scalar condition of (c) gives the true set as unique maximum.
  bediereffect    Prop. bediereffect: P(T is a maximiser) <= (1 - m_tau)^L by simulation.
  bedier          Thm bedier: KL bound >= simulated error of the thresholded estimator; second form <= first.
  rival           Prop. rival: conditional law of w(sigma) given the sum is binomial (randomised PIT, chi-square).
  dimfree         Prop. dimfree: display probability of a false split <= 2 max(delta, 2 delta^2), whatever n.
  dominance       Cor. dominance: second-order rooted support of each cluster <= unrooted support of the split.
  modelglobal     Cor. modelglobal: under its hypothesis, global separation of the exact irreversible-model supports.
  lba             Cor. lba (and Prop. quartet): formulas = enumeration; consistency iff (1-r)p1p3 < r(2-p1-p3), exactly.

Usage: audit_new_checks.py [OUT.json]   (default data/audit_new_checks.json; exits non-zero on any failure)"""
import sys, os, json, math, itertools, random, time
import numpy as np
import networkx as nx
from scipy.stats import chi2, binom


# ------------------------------------------------------------------ splits
def canon(A, X):
    A = frozenset(A); return A if 0 not in A else frozenset(X) - A


def compatible(A, B, X):
    Ac, Bc = X - A, X - B
    return not (A & B and A & Bc and Ac & B and Ac & Bc)


def all_splits(n):
    return [frozenset(A) for r in range(2, n - 1) for A in itertools.combinations(range(1, n), r)]


def compatible_sets(C, X):
    out = [frozenset()]
    def rec(i, cur):
        for j in range(i, len(C)):
            if all(compatible(C[j], t, X) for t in cur):
                nxt = cur + [C[j]]; out.append(frozenset(nxt)); rec(j + 1, nxt)
    rec(0, []); return out


class Space:
    """All non-trivial splits of n witnesses and all their compatible sets, as 0/1 matrices."""
    _cache = {}
    def __new__(cls, n):
        if n in cls._cache: return cls._cache[n]
        self = super().__new__(cls); self.n = n; self.X = frozenset(range(n))
        self.S = all_splits(n); self.ix = {s: i for i, s in enumerate(self.S)}
        self.sets = compatible_sets(self.S, self.X); self.set_ix = {s: i for i, s in enumerate(self.sets)}
        self.M = np.zeros((len(self.sets), len(self.S)))
        for i, C in enumerate(self.sets):
            for s in C: self.M[i, self.ix[s]] = 1
        self.size = self.M.sum(1)
        cls._cache[n] = self; return self

    def vec(self, d):
        return np.array([d.get(s, 0.0) for s in self.S])

    def gamma(self, m, T):
        """margin by the ratio definition, all compatible S != T"""
        T = frozenset(T); best = math.inf
        for S in self.sets:
            if S == T: continue
            A, B = T - S, S - T
            if not A: return -math.inf          # S strictly contains T (cannot happen for binary T)
            best = min(best, (sum(m[self.ix[s]] for s in A) - sum(m[self.ix[s]] for s in B)) / len(A))
        return best


def select(W, sp, t):
    """W: reps x splits supports. Returns (unique: S_hat = T for every selection rule, is_max: T is a maximiser).
    Only candidates (positive weight) may be selected."""
    pos = (W > 0).astype(float)
    valid = (pos @ sp.M.T) == sp.size                      # all splits of the set are candidates
    val = np.where(valid, W @ sp.M.T, -np.inf)
    vt = val[:, t].copy(); val[:, t] = -np.inf; other = val.max(1)
    okT = valid[:, t]
    return okT & (vt > other), okT & (vt >= other)


# ------------------------------------------------------------------ trees and the attraction model
def rooted_binary(n, R):
    nodes = [frozenset({i}) for i in range(n)]; cl = []
    while len(nodes) > 2:
        a, b = R.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; cl.append(c); nodes.append(c)
    return [frozenset({i}) for i in range(n)] + cl          # every edge once (root's two child edges included)


def contract(clades, n, k, R):
    """contract k random internal edges (clades of size 2..n-1) to make polytomies"""
    internal = [c for c in clades if 2 <= len(c) <= n - 1]
    drop = set(R.sample(internal, min(k, len(internal))))
    return [c for c in clades if c not in drop]


def true_splits(clades, n):
    X = frozenset(range(n)); return {canon(c, X) for c in clades if 2 <= len(c) <= n - 2}


def partition_law(n, clades, rates, al):
    """exact law of the reading partition of one unit (all 3^|E| innovation outcomes)."""
    order = sorted(range(len(clades)), key=lambda k: -len(clades[k]))       # higher edges first, lowest wins
    paths = [[k for k in order if leaf in clades[k]] for leaf in range(n)]
    law = {}
    for out in itertools.product((0, 1, 2), repeat=len(clades)):            # 0 none, 1 attractive, 2 own
        pr = 1.0
        for k, o in enumerate(out):
            pr *= (1 - rates[k]) if o == 0 else rates[k] * (al if o == 1 else 1 - al)
        if pr == 0: continue
        val = []
        for leaf in range(n):
            v = 'o'
            for k in paths[leaf]:
                if out[k]: v = 'a' if out[k] == 1 else k
            val.append(v)
        cls = {}
        for leaf, v in enumerate(val): cls.setdefault(v, set()).add(leaf)
        key = frozenset(frozenset(c) for c in cls.values())
        law[key] = law.get(key, 0.0) + pr
    return law


def proposals(part, n):
    X = frozenset(range(n)); return frozenset(canon(c, X) for c in part if 2 <= len(c) <= n - 2)


def types_from_law(law, sp):
    agg = {}
    for part, p in law.items():
        k = proposals(part, sp.n); agg[k] = agg.get(k, 0.0) + p
    keys = list(agg); probs = np.array([agg[k] for k in keys]); probs /= probs.sum()
    inc = np.zeros((len(keys), len(sp.S)))
    for i, k in enumerate(keys):
        for s in k: inc[i, sp.ix[s]] = 1
    return keys, probs, inc


def attraction_instance(n, R, binary=True, contract_k=0, rate=(0.02, 0.3), pendant=None):
    clades = rooted_binary(n, R)
    if not binary: clades = contract(clades, n, contract_k, R)
    rates = [R.uniform(*rate) if (pendant is None or len(c) > 1) else R.uniform(*pendant) for c in clades]
    al = R.random()
    law = partition_law(n, clades, rates, al)
    sp = Space(n); keys, probs, inc = types_from_law(law, sp)
    return dict(n=n, clades=clades, rates=rates, alpha=al, law=law, sp=sp, probs=probs, inc=inc,
                m=probs @ inc, T=frozenset(true_splits(clades, n)))


def simulate(inst, L, reps, rng):
    counts = rng.multinomial(L, inst['probs'], size=reps)
    return counts @ inst['inc']


def tex_e(x):
    a, b = f'{x:.1e}'.split('e'); return f'{a}\\times10^{{{int(b)}}}'


def num(x):
    return f'{x:,}'.replace(',', '{,}')


def se(p, reps): return math.sqrt(max(p, 1.0 / reps) * (1 - p) / reps)


# ================================================================== checks
def check_exclusive(seed=1):
    R = random.Random(seed); units = pairs = clpairs = 0; ok = True; cross_conflicts = 0
    for _ in range(400):                                     # collations
        n = R.randint(4, 9); X = frozenset(range(n)); prev = []
        for _u in range(R.randint(5, 30)):
            O = frozenset(x for x in X if R.random() < 0.85)
            if len(O) < 4: continue
            k = R.randint(2, len(O)); r = {x: R.randrange(k) for x in O}
            cls = {}
            for x, v in r.items(): cls.setdefault(v, set()).add(x)
            cls = [frozenset(c) for c in cls.values()]
            props = list({c if min(O) not in c else O - c for c in cls if 2 <= len(c) <= len(O) - 2})
            for a, b in itertools.combinations(props, 2):
                pairs += 1; ok &= compatible(a, b, O)
            ref = R.choice(list(set(r.values())))
            clusters = [c for c in cls if r[next(iter(c))] != ref and len(c) >= 2]
            for a, b in itertools.combinations(clusters, 2):
                clpairs += 1; ok &= not (a & b) or a <= b or b <= a
            if O == X:
                for a in props:
                    for b in prev: cross_conflicts += not compatible(a, b, X)
                prev += props
            units += 1
    ok &= cross_conflicts > 0                               # the test is not vacuous: different units do conflict
    return dict(passed=bool(ok), numbers=dict(collations=400, units=units, split_pairs=pairs, cluster_pairs=clpairs,
                                               cross_unit_conflicts=cross_conflicts),
                instances=f'400 random collations, $4\\le n\\le9$, multistate units with random illegibility ({num(units)} units)',
                result=f'all {num(pairs)} same-unit split pairs compatible and all {num(clpairs)} cluster pairs laminar; '
                       f'{num(cross_conflicts)} conflicting pairs across units (control)')


def check_nphard(seed=2):
    R = random.Random(seed); ok = True; graphs = 0; sizes_ok = True; maxV = 0
    for _ in range(250):
        nv = R.randint(2, 10); H = nx.gnp_random_graph(nv, R.random(), seed=R.randrange(10**9)); maxV = max(maxV, nv)
        V, E = list(H.nodes), [tuple(e) for e in H.edges]
        W = [('p', v, i) for v in V for i in (0, 1)] + [('q', e, a, b) for e in E for a in (0, 1) for b in (0, 1)] + [('z', 1), ('z', 2)]
        X = frozenset(W); sizes_ok &= len(X) == 4 * len(E) + 2 * len(V) + 2
        side = {}
        for v in V:
            s = {('p', v, 0), ('p', v, 1)}
            for e in E:
                if v in e:
                    k = e.index(v); s |= {('q', e, a, b) for a in (0, 1) for b in (0, 1) if (a, b)[k] == 1}
            side[v] = frozenset(s)
        w = {v: R.randint(1, 20) for v in V}
        conf = {(a, b): not compatible(side[a], side[b], X) for a in V for b in V if a != b}
        best_H, best_S, opt_H, opt_S = -1, -1, set(), set()
        for r in range(nv + 1):
            for U in itertools.combinations(V, r):
                val = sum(w[v] for v in U)
                if all(not H.has_edge(a, b) for a, b in itertools.combinations(U, 2)):            # independent in H
                    if val > best_H: best_H, opt_H = val, {frozenset(U)}
                    elif val == best_H: opt_H.add(frozenset(U))
                if all(not conf[a, b] for a, b in itertools.combinations(U, 2)):                  # compatible splits
                    if val > best_S: best_S, opt_S = val, {frozenset(U)}
                    elif val == best_S: opt_S.add(frozenset(U))
        ok &= best_H == best_S and opt_H == opt_S; graphs += 1
    ok &= sizes_ok
    return dict(passed=bool(ok), numbers=dict(graphs=graphs, max_vertices=maxV),
                instances=f'{graphs} random graphs, $2\\le|V|\\le{maxV}$, integer weights $1$--$20$',
                result='MWIS of $H$ (brute force) equals the maximum-weight compatible set of the constructed splits, '
                       'same value and same optimal sets, in all; $4|E|+2|V|+2$ witnesses')


def tree_with_polytomies(n, R, p_contract):
    G = nx.Graph(); G.add_edges_from([(0, 'h'), (1, 'h'), (2, 'h')]); k = 0
    for leaf in range(3, n):
        u, v = R.choice(sorted(G.edges, key=str)); w = f'i{k}'; k += 1
        G.remove_edge(u, v); G.add_edges_from([(u, w), (w, v), (w, leaf)])
    for e in [e for e in G.edges if not isinstance(e[0], int) and not isinstance(e[1], int)]:
        if R.random() < p_contract and G.has_edge(*e): G = nx.contracted_nodes(G, e[0], e[1], self_loops=False)
    return nx.Graph(G)


def splits_of_graph(G, n):
    X = frozenset(range(n)); S = set()
    for u, v in G.edges:
        H = G.copy(); H.remove_edge(u, v)
        side = frozenset(x for x in nx.node_connected_component(H, u) if isinstance(x, int))
        if 2 <= len(side) <= n - 2: S.add(canon(side, X))
    return S


def tree_from_splits(S, n):
    """root at witness 0: the sides avoiding 0 are a laminar family; hang each on its smallest strict superset."""
    Y = frozenset(range(1, n)); nodes = sorted(set(S) | {Y}, key=len)
    G = nx.Graph(); G.add_edge(0, ('c', Y))
    for c in nodes:
        if c == Y: continue
        par = min((d for d in nodes if c < d), key=len); G.add_edge(('c', c), ('c', par))
    for x in range(1, n):
        par = min((d for d in nodes if x in d), key=len); G.add_edge(x, ('c', par))
    return G


def check_buneman(seed=3):
    R = random.Random(seed); ok = True; trees = sets = 0; nonbinary = 0
    nm = lambda a, b: a.get('lab') == b.get('lab')
    for _ in range(300):
        n = R.randint(4, 10); G = tree_with_polytomies(n, R, R.random() * 0.6); X = frozenset(range(n))
        S = splits_of_graph(G, n); trees += 1; nonbinary += len(S) < n - 3
        ok &= all(compatible(a, b, X) for a, b in itertools.combinations(S, 2))          # "only if"
        H = tree_from_splits(S, n)
        ok &= splits_of_graph(H, n) == S
        ok &= all(d >= 3 for v, d in H.degree if not isinstance(v, int)) and all(H.degree(x) == 1 for x in range(n))
        for g in (G, H):
            for v in g.nodes: g.nodes[v]['lab'] = v if isinstance(v, int) else None
        ok &= nx.is_isomorphic(G, H, node_match=nm)                                      # unique up to isomorphism
    for _ in range(300):                                                                 # "if": random compatible sets
        n = R.randint(4, 10); X = frozenset(range(n)); C = all_splits(n); S = []
        p = R.random()
        for s in R.sample(C, len(C)):
            if R.random() < p and all(compatible(s, t, X) for t in S): S.append(s)
        H = tree_from_splits(set(S), n); sets += 1
        ok &= splits_of_graph(H, n) == set(S) and all(d >= 3 for v, d in H.degree if not isinstance(v, int))
    return dict(passed=bool(ok), numbers=dict(trees=trees, nonbinary_trees=nonbinary, compatible_sets=sets),
                instances=f'{trees} random stemmata ($4\\le n\\le10$, {nonbinary} non-binary) and {sets} random compatible sets',
                result='split sets pairwise compatible; split set $\\to$ tree $\\to$ split set exact, internal degrees $\\ge3$, '
                       'rebuilt tree isomorphic to the original as a leaf-labelled tree, in all')


def check_threshold(seed=4):
    R = random.Random(seed); ok = True; inst = 0; unrestricted_wrong = 0
    for _ in range(150):
        n = R.randint(5, 7); sp = Space(n)
        clades = contract(rooted_binary(n, R), n, R.randint(0, n - 3), R); T = true_splits(clades, n)
        if not T: continue
        t = R.randint(0, 8)
        w = np.array([R.randint(t + 1, t + 10) if s in T else R.randint(0, t) for s in sp.S], float)
        W = np.where(w > t, w, 0.0)[None, :]
        uniq, _ = select(W, sp, sp.set_ix[frozenset(T)]); ok &= bool(uniq[0])
        u2, _ = select(w[None, :], sp, sp.set_ix[frozenset(T)]); unrestricted_wrong += not u2[0]
        inst += 1
    return dict(passed=bool(ok), numbers=dict(instances=inst, unrestricted_optimum_not_T=unrestricted_wrong),
                instances=f'{inst} random stemmata, $5\\le n\\le7$ (binary and not), integer weights meeting the hypothesis',
                result=f'restricted optimum is $\\T$, uniquely, in all (all compatible sets enumerated); '
                       f'without the threshold the optimum differs from $\\T$ in {unrestricted_wrong}')


def pairwise_bound(inst, L):
    m, sp, T = inst['m'], inst['sp'], inst['T']
    F = [i for i, s in enumerate(sp.S) if s not in T and m[i] > 0]
    tot = sum((1 - (math.sqrt(m[sp.ix[s]]) - math.sqrt(m[f])) ** 2) ** L
              for f in F for s in T if not compatible(sp.S[f], s, sp.X))
    return tot + sum((1 - m[sp.ix[s]]) ** L for s in T)


def first_L(f, target, Lmax=10**7):
    L = 1
    while f(L) > target and L < Lmax: L = int(L * 1.05) + 1
    lo, hi = max(1, int(L / 1.05) - 1), L
    while lo < hi:
        mid = (lo + hi) // 2
        if f(mid) <= target: hi = mid
        else: lo = mid + 1
    return lo


def check_pairwise(seed=5):
    R = random.Random(seed); rng = np.random.default_rng(seed); rows = []; ok = True
    for n, want in ((5, 4), (6, 2)):
        found = 0
        for _ in range(200):
            if found == want: break
            inst = attraction_instance(n, R); m, sp, T = inst['m'], inst['sp'], inst['T']
            F = [i for i, s in enumerate(sp.S) if s not in T and m[i] > 0]
            hyp = all(m[sp.ix[s]] > 0 for s in T) and all(m[f] < m[sp.ix[s]] for f in F for s in T if not compatible(sp.S[f], s, sp.X))
            if not hyp: continue
            found += 1; t = sp.set_ix[T]
            for target in (0.5, 0.05):
                L = first_L(lambda L: pairwise_bound(inst, L), target); b = pairwise_bound(inst, L)
                reps = 40000 if n == 5 else 20000
                uniq, _ = select(simulate(inst, L, reps, rng), sp, t); err = 1 - uniq.mean()
                passed = err - 4 * se(err, reps) <= b; ok &= passed
                rows.append(dict(n=n, alpha=round(inst['alpha'], 3), K=len(F), L=L, bound=round(b, 4), sim_error=round(err, 5),
                                 reps=reps, passed=bool(passed)))
                print('  pairwise', rows[-1], flush=True)
    ok &= len(rows) == 12
    ratio = [r['sim_error'] / r['bound'] for r in rows]
    return dict(passed=bool(ok), numbers=dict(rows=rows, max_sim_over_bound=round(max(ratio), 4)),
                instances=f'{len(rows) // 2} attraction-model instances meeting the hypotheses ($n=5,6$), two lengths each '
                          f'(bound $\\approx0.5$ and $\\approx0.05$, $L={min(r["L"] for r in rows)}$--${max(r["L"] for r in rows)}$), '
                          f'$2\\times10^4$--$4\\times10^4$ simulated collations',
                result=f'bound $\\ge$ simulated $\\Pr(\\hat S_L\\neq\\T)$ (ties counted as errors) in all; largest ratio '
                       f'simulated/bound ${max(ratio):.3f}$')


def bernstein(m, T, sp, g, L):
    return sum(2 * math.exp(-L * g * g / (8 * m[i] + 4 * g / 3)) for i, s in enumerate(sp.S) if s in T or m[i] > 0)


def hoeffding(m, T, sp, g, L):
    K = sum(1 for i, s in enumerate(sp.S) if s not in T and m[i] > 0)
    return 2 * (sp.n - 3 + K) * math.exp(-L * g * g / 2)


def check_consistency(seed=6):
    R = random.Random(seed); rng = np.random.default_rng(seed); rows = []; ok = True; tries = 0
    while len(rows) < 8 and tries < 400:
        tries += 1; n = 5 if len(rows) < 6 else 6
        inst = attraction_instance(n, R, rate=(0.08, 0.35)); m, sp, T = inst['m'], inst['sp'], inst['T']
        g = sp.gamma(m, T)
        if g < 0.03: continue
        for name, f in (('Bernstein', bernstein), ('Hoeffding', hoeffding)):
            L = first_L(lambda L: f(m, T, sp, g, L), 0.5)
            reps = 20000
            uniq, _ = select(simulate(inst, L, reps, rng), sp, sp.set_ix[T]); err = 1 - uniq.mean()
            bB, bH = bernstein(m, T, sp, g, L), hoeffding(m, T, sp, g, L)
            passed = err - 4 * se(err, reps) <= min(bB, bH); ok &= passed
            rows.append(dict(n=n, gamma=round(g, 4), at=name, L=L, bernstein=round(bB, 4), hoeffding=round(bH, 4),
                             sim_error=round(err, 5), passed=bool(passed)))
            print('  consistency(a)', rows[-1], flush=True)
    ok &= len(rows) == 8
    # (b) inconsistency: designed law on n = 5, a false tau crossing sigma1 with m_tau > m_sigma1 and compatible
    # with sigma2, so S = {tau, sigma2} has m(S) > m(T); and an attraction-model long-branch instance.
    sp = Space(5); X = sp.X
    s1, s2 = frozenset({1, 2}), frozenset({3, 4})                       # T = 0|12|34 shape: 12|034 and 34|012
    tau = next(s for s in sp.S if not compatible(s, s1, X) and compatible(s, s2, X))
    keys = [frozenset({s1}), frozenset({s2}), frozenset({tau}), frozenset({s1, s2}), frozenset()]
    probs = np.array([0.12, 0.20, 0.22, 0.06, 0.40]); inc = np.zeros((5, len(sp.S)))
    for i, k in enumerate(keys):
        for s in k: inc[i, sp.ix[s]] = 1
    des = dict(probs=probs, inc=inc, sp=sp, T=frozenset({s1, s2}), m=probs @ inc)
    lb = None
    for _ in range(3000):
        cand = attraction_instance(5, R, rate=(0.0, 0.05), pendant=(0.2, 0.6))
        mT = sum(cand['m'][cand['sp'].ix[s]] for s in cand['T'])
        best = max(sum(cand['m'][cand['sp'].ix[s]] for s in S) for S in cand['sp'].sets)
        if best > mT + 0.02: lb = cand; break
    incons = []
    for tag, inst in (('designed', des), ('long branch', lb)):
        if inst is None: ok = False; continue
        sp_, T = inst['sp'], inst['T']; mT = sum(inst['m'][sp_.ix[s]] for s in T)
        gap = max(sum(inst['m'][sp_.ix[s]] for s in S) for S in sp_.sets) - mT
        fr = []
        for L in (25, 100, 400, 1600, 6400):
            uniq, _ = select(simulate(inst, L, 4000, rng), sp_, sp_.set_ix[T]); fr.append(round(float(uniq.mean()), 4))
        passed = fr[-1] <= 0.01 and fr[-1] < fr[0] + 1e-12; ok &= passed
        incons.append(dict(case=tag, gap=round(gap, 4), L=[25, 100, 400, 1600, 6400], freq_T=fr, passed=bool(passed)))
        print('  consistency(b)', incons[-1], flush=True)
    return dict(passed=bool(ok), numbers=dict(part_a=rows, part_b=incons),
                instances=f'(a) {len(rows) // 2} attraction-model instances with $\\gamma={min(r["gamma"] for r in rows):.3f}$--'
                          f'${max(r["gamma"] for r in rows):.3f}$ ($n=5,6$), at the length where each bound first reaches $0.5$ '
                          f'($L={num(min(r["L"] for r in rows))}$--${num(max(r["L"] for r in rows))}$), $2\\times10^4$ collations; '
                          f'(b) a designed law and a long-branch instance with $m(S)>m(\\T)$, $L=25$--$6{{,}}400$',
                result='(a) both bounds $\\ge$ simulated error in all (largest simulated error '
                       f'${max(r["sim_error"] for r in rows):.4f}$); (b) $\\Pr(\\hat S_L=\\T)$ falls from '
                       + ' and '.join(f'${c["freq_T"][0]:.2f}$ to ${c["freq_T"][-1]:.2f}$' for c in incons))


def check_samplesize(seed=7):
    R = random.Random(seed); rng = np.random.default_rng(seed); ok = True; inst_n = 0; worst = 0.0; sims = []
    for i in range(60):
        n = 5 if i % 3 else 6
        inst = attraction_instance(n, R, rate=(0.05, 0.35)); m, sp, T = inst['m'], inst['sp'], inst['T']
        g = sp.gamma(m, T)
        if g <= 0: continue
        K = sum(1 for j, s in enumerate(sp.S) if s not in T and m[j] > 0)
        mmax = max(m[j] for j, s in enumerate(sp.S) if s in T or m[j] > 0)
        for eps in (0.5, 0.1, 0.01, 0.001):
            L = math.ceil((8 * mmax + 4 * g / 3) / g ** 2 * math.log(2 * (n - 3 + K) / eps))
            b = bernstein(m, T, sp, g, L); ok &= b <= eps; worst = max(worst, b / eps); inst_n += 1
            if eps == 0.1 and len(sims) < 4 and L <= 200000:
                uniq, _ = select(simulate(inst, L, 5000, rng), sp, sp.set_ix[T]); err = 1 - uniq.mean()
                ok &= err - 4 * se(err, 5000) <= eps
                sims.append(dict(n=n, gamma=round(g, 4), L=L, bernstein=round(b, 4), sim_error=round(err, 5)))
                print('  samplesize', sims[-1], flush=True)
    ok &= len(sims) == 4
    return dict(passed=bool(ok), numbers=dict(cases=inst_n, max_bound_over_eps=round(worst, 4), simulated=sims),
                instances=f'{inst_n} (instance, $\\varepsilon$) cases: attraction-model instances with $\\gamma>0$ ($n=5,6$), '
                          f'$\\varepsilon\\in\\{{0.5,0.1,0.01,0.001\\}}$; 4 simulated at $\\varepsilon=0.1$',
                result=f'at the $L$ of the formula the Bernstein bound is $\\le\\varepsilon$ in all (largest ratio ${worst:.3f}$); '
                       + ('no error in the simulations at $\\varepsilon=0.1$' if max(s['sim_error'] for s in sims) == 0 else
                          f'largest simulated error at that $L$ ${max(s["sim_error"] for s in sims):.4f}$ ($\\varepsilon=0.1$)'))


def block_supports(inst, J, k, rho, reps, rng):
    """J i.i.d. blocks of k units: with probability rho the block repeats one unit k times, else k independent units."""
    Jc = rng.binomial(J, rho, size=reps)
    cc = rng.multinomial(Jc, inst['probs']); cf = rng.multinomial(k * (J - Jc), inst['probs'])
    return (k * cc + cf) @ inst['inc']


def check_blocks(seed=8):
    R = random.Random(seed); rng = np.random.default_rng(seed); rows = []; ok = True; tries = 0
    while len(rows) < 8 and tries < 400:
        tries += 1
        inst = attraction_instance(5, R, rate=(0.1, 0.4)); m, sp, T = inst['m'], inst['sp'], inst['T']
        g = sp.gamma(m, T)
        if g < 0.05: continue
        K = sum(1 for j, s in enumerate(sp.S) if s not in T and m[j] > 0)
        bnd = lambda J: 2 * (sp.n - 3 + K) * math.exp(-J * g * g / 2)
        for rho, k in ((1.0, 5), (0.5, 10)):
            J = first_L(bnd, 0.5); reps = 10000
            uniq, _ = select(block_supports(inst, J, k, rho, reps, rng), sp, sp.set_ix[T]); err = 1 - uniq.mean()
            passed = err - 4 * se(err, reps) <= bnd(J); ok &= passed
            rows.append(dict(gamma=round(g, 4), K=K, rho=rho, k=k, J=J, bound=round(bnd(J), 4), sim_error=round(err, 5), passed=bool(passed)))
            print('  blocks', rows[-1], flush=True)
    ok &= len(rows) == 8
    # control: blocks of k = 20 identical units, short J; the i.i.d. bound with L = Jk ignores the dependence
    ctrl = []; R2 = random.Random(seed + 100)
    while len(ctrl) < 12:
        inst = attraction_instance(5, R2, rate=(0.1, 0.4)); m, sp, T = inst['m'], inst['sp'], inst['T']
        g = sp.gamma(m, T)
        if g < 0.05: continue
        K = sum(1 for j, s in enumerate(sp.S) if s not in T and m[j] > 0)
        for J in (10, 20, 40, 80, 160, 320):
            reps = 10000
            uniq, _ = select(block_supports(inst, J, 20, 1.0, reps, rng), sp, sp.set_ix[T]); err = 1 - uniq.mean()
            bnd, naive = 2 * (sp.n - 3 + K) * math.exp(-J * g * g / 2), 2 * (sp.n - 3 + K) * math.exp(-20 * J * g * g / 2)
            passed = err - 4 * se(err, reps) <= bnd; ok &= passed
            ctrl.append(dict(gamma=round(g, 4), J=J, k=20, bound=round(bnd, 4), iid_bound_L_Jk=naive, sim_error=round(err, 5),
                             iid_bound_violated=bool(naive < err - 4 * se(err, reps)), passed=bool(passed)))
            print('  blocks control', ctrl[-1], flush=True)
    nviol = sum(c['iid_bound_violated'] for c in ctrl)
    return dict(passed=bool(ok), numbers=dict(rows=rows, control=ctrl, iid_bound_violations=nviol),
                instances=f'{len(rows) // 2} attraction-model instances ($n=5$, $\\gamma={min(r["gamma"] for r in rows):.3f}$--'
                          f'${max(r["gamma"] for r in rows):.3f}$), blocks of $k=5$ identical units or $k=10$ units '
                          f'repeated with probability $1/2$, $J={min(r["J"] for r in rows)}$--${max(r["J"] for r in rows)}$ blocks, $10^4$ collations',
                result=f'bound $\\ge$ simulated error in all (at bound $\\approx0.5$, simulated error '
                       f'${max(r["sim_error"] for r in rows):.4f}$ at most; with $k=20$ identical units and $J=10$--$320$, simulated error up to '
                       f'${max(c["sim_error"] for c in ctrl):.3f}$, still below the bound); control: the i.i.d.\\ bound with $L=Jk$ '
                       f'is exceeded by the simulated error in {nviol} of {len(ctrl)} short-$J$ cases')


def check_gammablocks(seed=9):
    R = random.Random(seed); rng = np.random.default_rng(seed); sp = Space(5); ok = True; rows = []
    binT = [S for S in sp.sets if len(S) == 2]                         # the 15 binary stemmata on 5 witnesses
    A = np.zeros((len(binT), len(sp.sets), len(sp.S))); B = np.zeros_like(A); size = np.ones((len(binT), len(sp.sets)))
    mask = np.ones((len(binT), len(sp.sets)), bool)
    for a, T in enumerate(binT):
        for b, S in enumerate(sp.sets):
            if S == T: mask[a, b] = False; continue
            for s in T - S: A[a, b, sp.ix[s]] = 1
            for s in S - T: B[a, b, sp.ix[s]] = 1
            size[a, b] = len(T - S)
    def glo(lo, hi):                                                    # reps x |binT|
        v = (np.einsum('abs,rs->rab', A, lo) - np.einsum('abs,rs->rab', B, hi)) / size[None]
        return np.where(mask[None], v, np.inf).min(2)
    N = len(sp.S); insts = []
    while len(insts) < 2:
        inst = attraction_instance(5, R, rate=(0.1, 0.4))
        if sp.gamma(inst['m'], inst['T']) > 0.05: insts.append(inst)
    for _ in range(3000):                                               # a long-branch instance: T is not the maximum
        cand = attraction_instance(5, R, rate=(0.0, 0.05), pendant=(0.2, 0.6))
        if sp.gamma(cand['m'], cand['T']) < -0.01: insts.append(cand); break
    for ii, inst in enumerate(insts):
        m = inst['m']; gam = glo(m[None], m[None])[0]                   # gamma of every binary stemma
        for delta in (0.05, 0.2):
            for J, k, rho in ((50, 4, 1.0), (400, 4, 1.0), (3200, 4, 0.5)):
                reps = 2000; L = J * k; h = math.sqrt(math.log(2 * N / delta) / (2 * J))
                W = block_supports(inst, J, k, rho, reps, rng); f = W / L
                gl = glo(f - h, f + h)
                cover = float((gl <= gam[None] + 1e-12).all(1).mean())
                pos = (W > 0).astype(float); valid = (pos @ sp.M.T) == sp.size
                val = np.where(valid, W @ sp.M.T, -np.inf); Sh = val.argmax(1)
                cert = np.zeros(reps, bool); false_cert = np.zeros(reps, bool); cert_true = np.zeros(reps, bool)
                for r in range(reps):
                    S = sp.sets[Sh[r]]
                    if len(S) != 2: continue
                    a = binT.index(S)
                    if gl[r, a] > 0:
                        cert[r] = True; false_cert[r] = gam[a] <= 0; cert_true[r] = S == inst['T']
                passed = cover >= 1 - delta - 3 * se(1 - delta, reps) and false_cert.mean() <= delta + 3 * se(delta, reps)
                ok &= passed
                rows.append(dict(instance=ii, gamma_T=round(float(gam[binT.index(inst['T'])]), 4), delta=delta, J=J, k=k, rho=rho,
                                 coverage=cover, certified=float(cert.mean()), false_certificates=float(false_cert.mean()),
                                 passed=bool(passed)))
                print('  gammablocks', rows[-1], flush=True)
    ok &= len(insts) == 3
    return dict(passed=bool(ok), numbers=dict(rows=rows),
                instances=f'3 attraction-model instances on $n=5$ (two with $\\gamma>0$, one long-branch with $\\gamma<0$), '
                          f'$\\delta=0.05,0.2$, $J=50,400,3{{,}}200$ blocks of $4$ units (identical, or repeated with probability $1/2$), '
                          f'$2{{,}}000$ collations each',
                result=f'$\\gamma_{{\\mathrm{{lo}}}}(T)\\le\\gamma(T)$ for all $15$ binary stemmata simultaneously '
                       + ('in every simulated collation' if min(r['coverage'] for r in rows) == 1 else
                          f'in a fraction ${min(r["coverage"] for r in rows):.3f}$ of collations at least')
                       + ' (required $1-\\delta$); '
                       + ('no false certificate' if max(r['false_certificates'] for r in rows) == 0 else
                          f'false certificates ${max(r["false_certificates"] for r in rows):.3f}$ at most (allowed $\\delta$)')
                       + f'; the true stemma is certified in a fraction ${max(r["certified"] for r in rows if r["J"] == 50):.2f}$ '
                       f'or less at $J=50$ and up to ${max(r["certified"] for r in rows):.2f}$ at $J=3{{,}}200$')


def illegibility_law(n, delta, R, style):
    """a joint law of the illegible set I with every marginal <= delta (max marginal = delta)"""
    subs = [frozenset(c) for r in range(1, n + 1) for c in itertools.combinations(range(n), r)]
    if style == 'singleton':
        q = {frozenset({R.randrange(n)}): 1.0}
    elif style == 'block':
        q = {frozenset(R.sample(range(n), R.randint(2, n))): 1.0}
    else:
        chosen = R.sample(subs, R.randint(1, 6)); g = np.random.default_rng(R.randrange(10**9)).dirichlet([0.5] * len(chosen))
        q = dict(zip(chosen, g))
    marg = [sum(p for I, p in q.items() if i in I) for i in range(n)]
    c = delta / max(marg); q = {I: c * p for I, p in q.items()}; q[frozenset()] = 1 - sum(q.values())
    assert q[frozenset()] >= -1e-12
    return q


def restricted_supports(law, n, delta, R, laws=None):
    """exact m and mt when the illegibility law depends on the readings (one law per reading partition)."""
    X = frozenset(range(n)); S = all_splits(n); m = {s: 0.0 for s in S}; mt = {s: 0.0 for s in S}
    laws = laws if laws is not None else {}
    for part, pr in law.items():
        for s in proposals(part, n): m[s] += pr
        if part not in laws: laws[part] = (R.choice(('mix', 'mix', 'singleton', 'block')), R.randrange(10**9))
        style, sd = laws[part]; q = illegibility_law(n, delta, random.Random(sd), style)
        for I, pi in q.items():
            if pi <= 0: continue
            O = X - I
            pp = {(c & O, O - c) for c in part if len(c & O) >= 2 and len(O - c) >= 2}
            for s in S:
                a, b = s & O, (X - s) & O
                if len(a) >= 2 and len(b) >= 2 and ((a, b) in pp or (b, a) in pp): mt[s] += pr * pi
    return m, mt, laws


def targeted_supports(law, n, delta, sprop, I):
    """illegibility that depends on the readings: when the unit proposes sprop, the witnesses I are illegible together
    with probability delta (each witness's conditional rate is delta); otherwise all are legible."""
    X = frozenset(range(n)); S = all_splits(n); m = {s: 0.0 for s in S}; mt = {s: 0.0 for s in S}
    for part, pr in law.items():
        props = proposals(part, n)
        for s in props: m[s] += pr
        q = {I: delta, frozenset(): 1 - delta} if sprop in props else {frozenset(): 1.0}
        for II, pi in q.items():
            O = X - II
            pp = {(c & O, O - c) for c in part if len(c & O) >= 2 and len(O - c) >= 2}
            for s in S:
                a, b = s & O, (X - s) & O
                if len(a) >= 2 and len(b) >= 2 and ((a, b) in pp or (b, a) in pp): mt[s] += pr * pi
    return m, mt


def check_textillegible(seed=10):
    R = random.Random(seed); ok = True; inst = 0; sharp_viol = 0; forced = 0; min_lo = min_hi = math.inf
    for trial in range(36):
        n = 5 if trial < 30 else 6; X = frozenset(range(n)); clades = rooted_binary(n, R)
        rates = [R.uniform(0.02, 0.25) for _ in clades]; al = R.random(); law = partition_law(n, clades, rates, al)
        delta = R.uniform(0.01, 1.0 / n); T = true_splits(clades, n); S = all_splits(n)
        m, mt, laws = restricted_supports(law, n, delta, R)
        mu = sum(m.values()); d = lambda s, t: min(len(s ^ t), len(s ^ (X - t)))
        for s in S:
            lo = mt[s] - (1 - 4 * delta) * m[s]; hi = m[s] + 2 * delta * mu - mt[s]
            min_lo, min_hi = min(min_lo, lo), min(min_hi, hi); ok &= lo >= -1e-12 and hi >= -1e-12
            sharp_viol += mt[s] > m[s] + 2 * sum(m[x] * delta ** d(x, s) for x in S if x != s) + 1e-12
        inst += 1
        conflicts = [(t, s) for t in S if t not in T for s in T if not compatible(t, s, X)]
        rho = min(m[s] - m[t] for t, s in conflicts)
        if rho > 0:                                     # (c): delta just inside rho/(4+2mu), same law styles
            df = 0.9 * rho / (4 + 2 * mu); _, mtf, _ = restricted_supports(law, n, df, R, laws)
            vals = sorted(((sum(mtf[x] for x in C), C) for C in Space(n).sets), key=lambda z: -z[0])
            u = vals[0][1] == frozenset(T) and vals[1][0] < vals[0][0]; ok &= u; forced += 1
    # targeted laws on n = 6: illegible together exactly the witnesses that turn a proposed true split into a false one
    tgt = 0
    for trial in range(12):
        n = 6; X = frozenset(range(n)); clades = rooted_binary(n, R); S = all_splits(n); T = true_splits(clades, n)
        rates = [R.uniform(0.02, 0.25) if len(c) > 1 else R.uniform(0.0, 0.03) for c in clades]; al = R.random()
        law = partition_law(n, clades, rates, al); delta = R.uniform(0.05, 0.2)
        d = lambda s, t: min(len(s ^ t), len(s ^ (X - t)))
        pairs = [(s, t) for s in T for t in S if t not in T and d(s, t) == 2]
        s0, t0 = R.choice(pairs); I = frozenset(s0 ^ t0 if len(s0 ^ t0) == 2 else s0 ^ (X - t0))
        m, mt = targeted_supports(law, n, delta, s0, I); mu = sum(m.values())
        for s in S:
            lo = mt[s] - (1 - 4 * delta) * m[s]; hi = m[s] + 2 * delta * mu - mt[s]
            min_lo, min_hi = min(min_lo, lo), min(min_hi, hi); ok &= lo >= -1e-12 and hi >= -1e-12
            v = mt[s] > m[s] + 2 * sum(m[x] * delta ** d(x, s) for x in S if x != s) + 1e-12
            sharp_viol += v; tgt += v
        inst += 1
    return dict(passed=bool(ok), numbers=dict(instances=inst, min_slack_lower=min_lo, min_slack_upper=min_hi,
                                               sharp_bound_violations=sharp_viol, sharp_violations_targeted=tgt,
                                               condition_c_instances=forced),
                instances=f'{inst} random rooted trees ($n=5$: 30, $n=6$: 18), attraction model, all innovation outcomes; for each '
                          'reading partition its own joint law of the illegible set (singletons, blocks, random mixtures, dependent across '
                          'witnesses; on 12 trees, two witnesses illegible together exactly when they turn a proposed true split into a '
                          'false one), every conditional rate $\\le\\delta$',
                result=f'$\\tilde m_\\sigma\\ge(1-4\\delta)m_\\sigma$ and $\\tilde m_\\tau\\le m_\\tau+2\\delta\\mu$ for every split '
                       f'(smallest slacks ${tex_e(min_lo)}$, ${tex_e(min_hi)}$); the independent-case bound with $\\delta^{{d}}$ fails for '
                       f'{sharp_viol} splits (control; {tgt} of them under the targeted laws); in all {forced} instances with $\\varrho>0$, rerun at $\\delta=0.9\\varrho/(4+2\\mu)$, '
                       'the true set is the unique maximum of $\\tilde m$')


def polytomy_instance(R, n):
    for _ in range(500):
        inst = attraction_instance(n, R, binary=False, contract_k=R.randint(1, n - 3), rate=(0.03, 0.3))
        T, sp, m = inst['T'], inst['sp'], inst['m']
        if not T or len(T) == n - 3: continue
        comp = [i for i, s in enumerate(sp.S) if s not in T and m[i] > 0 and all(compatible(s, t, sp.X) for t in T)]
        if comp: return inst, comp
    raise RuntimeError('no polytomy instance')


def check_bediereffect(seed=11):
    R = random.Random(seed); rng = np.random.default_rng(seed); rows = []; ok = True
    for i in range(6):
        inst, comp = polytomy_instance(R, 5 if i < 4 else 6); sp, m = inst['sp'], inst['m']
        mt = max(m[j] for j in comp); t = sp.set_ix[inst['T']]
        for target in (0.5, 0.1):
            L = max(1, math.ceil(math.log(target) / math.log(1 - mt))); b = (1 - mt) ** L; reps = 20000
            uniq, ismax = select(simulate(inst, L, reps, rng), sp, t); p = float(ismax.mean())
            passed = p - 4 * se(p, reps) <= b; ok &= passed
            rows.append(dict(n=sp.n, splits_T=len(inst['T']), m_tau=round(mt, 4), L=L, bound=round(b, 4),
                             P_T_is_a_maximiser=round(p, 5), P_unique=round(float(uniq.mean()), 5), passed=bool(passed)))
            print('  bediereffect', rows[-1], flush=True)
    return dict(passed=bool(ok), numbers=dict(rows=rows),
                instances='6 attraction-model trees with a polytomy ($n=5,6$), the largest $m_\\tau$ among false splits compatible with '
                          f'$\\T$ (${min(r["m_tau"] for r in rows):.3f}$--${max(r["m_tau"] for r in rows):.3f}$), two lengths each '
                          f'($L={min(r["L"] for r in rows)}$--${max(r["L"] for r in rows)}$), $2\\times10^4$ collations',
                result='$\\Pr(\\T$ is a maximiser$)\\le(1-m_\\tau)^L$ in all (ratio simulated/bound '
                       f'${min(r["P_T_is_a_maximiser"] / r["bound"] for r in rows):.2f}$--${max(r["P_T_is_a_maximiser"] / r["bound"] for r in rows):.2f}$)')


def kl(x, y):
    t = 0.0
    if x > 0: t += x * math.log(x / y)
    if x < 1: t += (1 - x) * math.log((1 - x) / (1 - y))
    return t


def check_bedier(seed=12):
    R = random.Random(seed); rng = np.random.default_rng(seed); rows = []; ok = True; second_ok = True; tries = 0; found = 0
    while found < 5 and tries < 300:
        tries += 1
        inst, _ = polytomy_instance(R, 5 if found < 4 else 6); sp, m, T = inst['sp'], inst['m'], inst['T']
        b = min(m[sp.ix[s]] for s in T); a = max(m[i] for i, s in enumerate(sp.S) if s not in T)
        if not b > a: continue
        found += 1; F = [i for i, s in enumerate(sp.S) if s not in T and m[i] > 0]
        bound = lambda L, c: sum(math.exp(-L * kl(c, m[sp.ix[s]])) for s in T) + sum(math.exp(-L * kl(c, m[i])) for i in F)
        cstar = ((math.sqrt(a) + math.sqrt(b)) / 2) ** 2
        for L in (10, 30, 100, 300, 1000, 3000, 10000):                # second form: at c*, bound <= (|T|+K) e^{-L(sqrt b - sqrt a)^2/4}
            second_ok &= bound(L, cstar) <= (len(T) + len(F)) * math.exp(-L * (math.sqrt(b) - math.sqrt(a)) ** 2 / 4) + 1e-12
        for c in (cstar, a + R.uniform(0.1, 0.9) * (b - a)):
            L = first_L(lambda L: bound(L, c), 0.5); reps = 20000
            W = simulate(inst, L, reps, rng); W = np.where(W > c * L, W, 0.0)
            uniq, _ = select(W, sp, sp.set_ix[T]); err = 1 - float(uniq.mean())
            passed = err - 4 * se(err, reps) <= bound(L, c); ok &= passed
            rows.append(dict(n=sp.n, a=round(a, 4), b=round(b, 4), c=round(c, 4), L=L, bound=round(bound(L, c), 4),
                             sim_error=round(err, 5), passed=bool(passed)))
            print('  bedier', rows[-1], flush=True)
    ok &= second_ok and found == 5
    return dict(passed=bool(ok), numbers=dict(rows=rows, second_form_holds=bool(second_ok)),
                instances=f'5 attraction-model trees with a polytomy ($n=5,6$) and $b>a$, $c=c^\\ast$ and a random $c\\in(a,b)$, '
                          f'at the length where the bound first reaches $0.5$ ($L={num(min(r["L"] for r in rows))}$--${num(max(r["L"] for r in rows))}$), '
                          '$2\\times10^4$ collations',
                result=f'KL bound $\\ge$ simulated $\\Pr(\\hat S_{{L,c}}\\neq\\T)$ in all (largest simulated error '
                       f'${max(r["sim_error"] for r in rows):.4f}$); at $c^\\ast$ the KL bound is below the second form at every '
                       '$L$ tested ($10$--$10^4$)')


def check_rival(seed=13):
    R = random.Random(seed); rng = np.random.default_rng(seed); rows = []; ok = True; both = 0.0
    for i in range(4):
        inst = attraction_instance(5 if i < 3 else 6, R, rate=(0.05, 0.35)); sp, m, T = inst['sp'], inst['m'], inst['T']
        s = R.choice(sorted(T, key=sorted)); si = sp.ix[s]
        rivals = [j for j, t in enumerate(sp.S) if not compatible(t, s, sp.X) and m[j] > 0]
        if not rivals: continue
        ti = max(rivals, key=lambda j: m[j])
        both = max(both, float((inst['probs'] * inst['inc'][:, si] * inst['inc'][:, ti]).sum()))   # exact: never both
        p = m[si] / (m[si] + m[ti]); L = 200; reps = 40000
        W = simulate(inst, L, reps, rng); k, Nn = W[:, si], W[:, si] + W[:, ti]
        def pit_chi2(pp):
            V = rng.random(reps); u = np.where(Nn > 0, binom.cdf(k - 1, Nn, pp) + V * binom.pmf(k, Nn, pp), V)
            h = np.histogram(u, bins=10, range=(0, 1))[0]; e = reps / 10
            return float(chi2.sf(((h - e) ** 2 / e).sum(), 9))
        p_ok, p_wrong = pit_chi2(p), pit_chi2(min(0.99, p * 1.15))
        passed = p_ok > 1e-3 and p_wrong < 1e-6; ok &= passed
        rows.append(dict(n=sp.n, m_sigma=round(m[si], 4), m_tau=round(m[ti], 4), L=L, reps=reps, p_value=p_ok,
                         p_value_wrong_law=p_wrong, passed=bool(passed)))
        print('  rival', rows[-1], flush=True)
    ok &= both == 0.0 and len(rows) >= 3
    return dict(passed=bool(ok), numbers=dict(rows=rows, max_prob_unit_proposes_both=both),
                instances=f'{len(rows)} attraction-model trees ($n=5,6$), a true split and its strongest crossing rival, '
                          '$L=200$, $4\\times10^4$ collations',
                result='exact unit law: no unit proposes both; randomised probability-integral transform of $w(\\sigma)$ given the sum, '
                       f'chi-square on deciles: $p={min(r["p_value"] for r in rows):.2f}$--${max(r["p_value"] for r in rows):.2f}$ under the '
                       f'binomial law, $p<10^{{-6}}$ for the control with the ratio raised by $15\\%$')


def check_dimfree(seed=14):
    """prop:dimfree: for binary T' and sigma'=A|B not in Sigma(T'), with independent illegibility delta,
    P(sigma'|O non-trivial and displayed by T') <= 2 max(delta, 2 delta^2); and, structurally, every displayed
    pattern leaves one of the four parts A1, A2, B1, B2 of the proof wholly illegible (which gives 4 delta under
    the conditional hypothesis)."""
    import networkx as nx
    R = random.Random(seed); ok = True; worst = 0.0; trees = 0; splits_tested = 0; struct_fail = 0
    for n in (5, 6, 7, 8):
        X = frozenset(range(n))
        for trial in range(20 if n < 8 else 8):
            T = nx.Graph(); T.add_edges_from([(0, 'r'), (1, 'r'), (2, 'r')]); k = 0
            for leaf in range(3, n):
                u, v = R.choice(list(T.edges())); mid = ('i', k); k += 1
                T.remove_edge(u, v); T.add_edges_from([(u, mid), (mid, v), (mid, leaf)])
            trees += 1
            sides = []
            for u, v in T.edges():
                H = T.copy(); H.remove_edge(u, v)
                c = frozenset(x for x in nx.node_connected_component(H, u) if isinstance(x, int))
                if 2 <= len(c) <= n - 2: sides += [c, X - c]
            def displayed(a, b):
                return any((a <= c and b <= X - c) for c in sides)
            def hull(S):
                S = list(S); E = set()
                for i in range(len(S)):
                    for j in range(i + 1, len(S)): E |= set(nx.shortest_path(T, S[i], S[j]))
                return E if len(S) > 1 else set(S)
            for r in range(2, n - 1):
                for A in itertools.combinations(range(n), r):
                    A = frozenset(A); B = X - A
                    if 0 not in A or any(A == c for c in sides): continue
                    splits_tested += 1
                    I = hull(A) & hull(B)
                    # leaf u of I (internal vertex) and the four parts
                    Isub = T.subgraph(I); u = next(w for w in I if Isub.degree(w) <= 1)
                    parts = []
                    for nb in T.neighbors(u):
                        H = T.copy(); H.remove_node(u)
                        comp = frozenset(x for x in nx.node_connected_component(H, nb) if isinstance(x, int)); parts.append(comp)
                    pats = []
                    for k2 in range(n + 1):
                        for ill in itertools.combinations(range(n), k2):
                            O = X - frozenset(ill); a, b = A & O, B & O
                            if len(a) >= 2 and len(b) >= 2 and displayed(a, b):
                                pats.append(k2)
                                Ai = [p & A for p in parts if p & A]; Bi = [p & B for p in parts if p & B]
                                if not any(q <= frozenset(ill) for q in Ai + Bi): struct_fail += 1
                    for d in (0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9):
                        prob = sum(d ** k2 * (1 - d) ** (n - k2) for k2 in pats)
                        bound = 2 * max(d, 2 * d * d); worst = max(worst, prob / bound)
                        ok &= prob <= bound + 1e-12
    ok &= struct_fail == 0
    return dict(passed=bool(ok), numbers=dict(trees=trees, false_splits=splits_tested, worst_ratio=worst, structural_failures=struct_fail),
                instances=f'{trees} random binary trees $T\'$ ($n=5$--$8$), all {splits_tested} false splits, all legibility patterns, '
                          '$\\delta\\in\\{0.01,\\dots,0.9\\}$',
                result=f'display probability $\\le2\\max(\\delta,2\\delta^2)$ in all (largest ratio {worst:.3f}); every displayed pattern '
                       'leaves one of the four parts of the proof wholly illegible')


# ------------------------------------------------------------------ irreversible model (alpha = 1), Section 8
def irreversible_m(clades, rates, n):
    """exact supports of the irreversible model: all 2^|E| mark patterns; Y = union of marked clades proposes Y|rest"""
    full = (1 << n) - 1; masks = [sum(1 << i for i in c) for c in clades]; m = {}
    for pat in range(1 << len(clades)):
        pr = 1.0; Y = 0
        for k, q in enumerate(rates):
            if pat >> k & 1: pr *= q; Y |= masks[k]
            else: pr *= 1 - q
        c = bin(Y).count('1')
        if 2 <= c <= n - 2:
            if Y & 1: Y = full ^ Y
            s = frozenset(i for i in range(n) if Y >> i & 1); m[s] = m.get(s, 0.0) + pr
    return m


def check_modelglobal(seed=15):
    """cor:modelglobal: binary T, rates meeting min_e pi_e prod_f (1-pi_f) > sum_{disjoint pairs} pi_e pi_f;
    exact m by enumeration of all mark patterns: min true m > max false m (and gamma > 0 for n <= 6). Rooted on an
    edge (two halves) or at an internal vertex. Instances that miss the hypothesis are also evaluated, to show the
    separation test can fail."""
    R = random.Random(seed); ok = True; acc = rej = rej_fail = 0; worst = math.inf; gam_min = math.inf; tried = 0
    target = {5: 60, 6: 60, 7: 40, 8: 25}
    for n, N in target.items():
        X = frozenset(range(n)); got = 0
        while got < N:
            tried += 1
            clades = rooted_binary(n, R)
            roots = [c for c in clades if (X - c) in clades and len(c) >= 2]
            if R.random() < 0.5 and roots:
                drop = R.choice(roots); clades = [c for c in clades if c != drop]          # root at a vertex
            lo = R.choice((1e-4, 1e-3, 1e-2))
            rates = [math.exp(R.uniform(math.log(0.01), math.log(0.6))) if len(c) > 1
                     else math.exp(R.uniform(math.log(lo), math.log(0.2))) for c in clades]
            E = len(clades); keep = math.prod(1 - q for q in rates)
            lhs = min(rates[k] * keep for k in range(E) if 2 <= len(clades[k]) <= n - 2)
            rhs = sum(rates[e] * rates[f] for e, f in itertools.combinations(range(E), 2) if not clades[e] & clades[f])
            hyp = lhs > rhs
            if not hyp and rej >= 4 * acc + 20: continue          # keep a bounded sample of the rest
            m = irreversible_m(clades, rates, n); T = true_splits(clades, n)
            tmin = min(m.get(s, 0.0) for s in T); fmax = max([m.get(s, 0.0) for s in all_splits(n) if canon(s, X) not in T] + [0.0])
            sep = tmin > fmax
            if hyp:
                acc += 1; got += 1; ok &= sep; worst = min(worst, tmin - fmax)
                ok &= tmin >= lhs - 1e-15 and sum(v for s, v in m.items() if s not in T) <= rhs + 1e-15
                if n <= 6:
                    sp = Space(n); g = sp.gamma(sp.vec(m), T); gam_min = min(gam_min, g); ok &= g > 0
            else:
                rej += 1; rej_fail += not sep
    ok &= rej_fail > 0                                       # the test is able to fail outside the hypothesis
    return dict(passed=bool(ok), numbers=dict(meeting=acc, not_meeting=rej, not_meeting_separation_fails=rej_fail,
                                             smallest_margin=worst, smallest_gamma_n_le_6=gam_min),
                instances=f'{acc} random rooted binary trees ($n=5$--$8$, rooted on an edge or at a vertex) with rates meeting the '
                          'hypothesis; all $2^{|\\mathcal E|}$ mark patterns enumerated',
                result=f'global separation in all (smallest $\\min_\\T m-\\max_\\F m={tex_e(worst)}$), $\\gamma>0$ for $n\\le6$; '
                       f'of {rej} instances missing the hypothesis, separation fails in {rej_fail}')


def quartet_enum(p, r1, r2, one=1):
    """irreversible quartet T=12|34, all 2^6 mark patterns, exact in the number type of the inputs"""
    clades = [{0}, {1}, {2}, {3}, {0, 1}, {2, 3}]; rates = list(p) + [r1, r2]; m = {12: 0 * one, 13: 0 * one, 14: 0 * one}
    key = {frozenset({0, 1}): 12, frozenset({2, 3}): 12, frozenset({0, 2}): 13, frozenset({1, 3}): 13,
           frozenset({0, 3}): 14, frozenset({1, 2}): 14}
    for pat in itertools.product((0, 1), repeat=6):
        pr = one; Y = set()
        for c, q, b in zip(clades, rates, pat):
            pr *= q if b else one - q
            if b: Y |= c
        if len(Y) == 2: m[key[frozenset(Y)]] += pr
    return m


def quartet_formulas(p, r1, r2, one=1):
    b = lambda x: one - x; p1, p2, p3, p4 = p
    return {12: b(r2) * b(p3) * b(p4) * (r1 + b(r1) * p1 * p2) + b(r1) * b(p1) * b(p2) * (r2 + b(r2) * p3 * p4),
            13: b(r1) * b(r2) * (p1 * p3 * b(p2) * b(p4) + p2 * p4 * b(p1) * b(p3)),
            14: b(r1) * b(r2) * (p1 * p4 * b(p2) * b(p3) + p2 * p3 * b(p1) * b(p4))}


def check_lba(seed=16):
    """cor:lba (with prop:quartet): exact rational arithmetic. (1) the formulas of prop:quartet equal enumeration of the
    2^6 mark patterns, at random general rates and at every point below; (2) with p2=p4=0, r1=r2=r<1: the reduced
    supports, and m12 > max(m13, m14) iff (1-r)p1p3 < r(2-p1-p3), m13 > max(m12, m14) iff the reverse strict inequality,
    on a grid, at random points, and at points within 1e-9 of the boundary on either side."""
    from fractions import Fraction as Fr
    R = random.Random(seed); one = Fr(1); ok = True; formula_pts = 0
    rnd = lambda: Fr(R.randint(0, 10 ** 6), 10 ** 6)
    for _ in range(300):                                     # general rates
        p = [rnd() for _ in range(4)]; r1, r2 = rnd(), rnd()
        ok &= quartet_enum(p, r1, r2, one) == quartet_formulas(p, r1, r2, one); formula_pts += 1
    pts = [(Fr(i, 10), Fr(j, 10), Fr(k, 10)) for i in range(10) for j in range(11) for k in range(11)]   # r, p1, p3
    pts += [(Fr(R.randint(0, 10 ** 6 - 1), 10 ** 6), rnd(), rnd()) for _ in range(1500)]
    near = 0
    while near < 1000:
        r = Fr(R.randint(1, 10 ** 6 - 1), 10 ** 6); p1 = rnd()
        p3 = r * (2 - p1) / ((1 - r) * p1 + r)               # boundary
        if not 0 < p3 < 1: continue
        for d in (Fr(1, 10 ** 9), -Fr(1, 10 ** 9)):
            if 0 <= p3 + d <= 1: pts.append((r, p1, p3 + d)); near += 1
        pts.append((r, p1, p3))
    side = dict(consistent=0, lba=0, boundary=0); mismatch = 0
    for r, p1, p3 in pts:
        p = (p1, 0 * one, p3, 0 * one)
        me = quartet_enum(p, r, r, one); mf = quartet_formulas(p, r, r, one); formula_pts += 1
        red = {12: r * (1 - r) * (2 - p1 - p3), 13: (1 - r) ** 2 * p1 * p3, 14: 0 * one}
        ok &= me == mf == red
        lhs, rhs = (1 - r) * p1 * p3, r * (2 - p1 - p3)
        cons = me[12] > max(me[13], me[14]); lba = me[13] > max(me[12], me[14])
        good = cons == (lhs < rhs) and lba == (lhs > rhs)
        mismatch += not good; ok &= good
        side['consistent' if lhs < rhs else 'lba' if lhs > rhs else 'boundary'] += 1
    ok &= side['consistent'] > 0 and side['lba'] > 0
    return dict(passed=bool(ok), numbers=dict(points=len(pts), near_boundary=near, formula_checks=formula_pts, mismatches=mismatch, **side),
                instances=f'{len(pts)} points $(r,p_1,p_3)$, $p_2=p_4=0$: a grid, random points and {near} within $10^{{-9}}$ of the '
                          f'boundary; and $300$ random general rate vectors; exact rational arithmetic',
                result=f'the formulas of Proposition~\\ref{{prop:quartet}} equal enumeration of the $2^6$ mark patterns at all '
                       f'{formula_pts} points; $m_{{12|34}}>\\max(m_{{13|24}},m_{{14|23}})$ exactly where $(1-r)p_1p_3<r(2-p_1-p_3)$ '
                       f'({side["consistent"]} points) and $13|24$ strictly heaviest exactly on the reverse side ({side["lba"]}); '
                       f'neither at the {side["boundary"]} points of equality; {mismatch} mismatches')


def second_order(n, clades, rates, al):
    """outcomes with at most two innovations of the attraction model: (prob, unrooted splits, rooted clusters)"""
    X = frozenset(range(n)); E = len(clades); order = sorted(range(E), key=lambda k: -len(clades[k])); out = []
    for kset in [()] + [(k,) for k in range(E)] + list(itertools.combinations(range(E), 2)):
        base = math.prod(rates[k] if k in kset else 1 - rates[k] for k in range(E))
        for rd in itertools.product('af', repeat=len(kset)):
            w = math.prod(al if v == 'a' else 1 - al for v in rd); read = dict(zip(kset, rd)); cls = {}
            for leaf in range(n):
                v = 'o'
                for k in order:
                    if leaf in clades[k] and k in read: v = 'a' if read[k] == 'a' else k
                cls.setdefault(v, set()).add(leaf)
            U = {canon(c, X) for c in map(frozenset, cls.values()) if 2 <= len(c) <= n - 2}
            Rt = {frozenset(c) for v, c in cls.items() if v != 'o' and 2 <= len(c) <= n - 1}
            out.append((base * w, U, Rt, kset))
    return out


def check_dominance(seed=17):
    """cor:dominance: random rooted trees (binary and with polytomies), attraction model, every outcome with at most two
    innovations enumerated. For every split tau=C|D not of the tree: Lambda^root(C), Lambda^root(D) <= Lambda^unr(tau);
    per outcome, a rooted proposal of C or D implies the unrooted proposal of tau; and the decomposition of the proof
    (disjoint pairs with union C: factor alpha^2 vs 1; nested pairs with C_e minus C_f = C: 1-alpha^2 in both; disjoint
    pairs with union D: unrooted only), computed from the clades alone, equals the enumerated second-order supports."""
    R = random.Random(seed); ok = True; trees = taus = strict = 0; worst = math.inf; outcome_fail = decomp_err = 0.0
    for trial in range(60):
        n = R.choice((5, 6, 7)); X = frozenset(range(n))
        clades = rooted_binary(n, R)
        if trial % 3 == 2: clades = contract(clades, n, R.randint(1, 2), R)
        rates = [R.uniform(0.02, 0.3) for _ in clades]; al = (0.0, 1.0, R.random(), R.random())[trial % 4]
        E = len(clades); outs = second_order(n, clades, rates, al); T = true_splits(clades, n); trees += 1
        LU, LR = {}, {}
        for pr, U, Rt, _ in outs:
            for s in U: LU[s] = LU.get(s, 0.0) + pr
            for c in Rt: LR[c] = LR.get(c, 0.0) + pr
            for c in Rt:
                if len(c) <= n - 2 and canon(c, X) not in T and canon(c, X) not in U: outcome_fail += 1
        def pair_w(e, f): return rates[e] * rates[f] * math.prod(1 - rates[g] for g in range(E) if g not in (e, f))
        dis, nest = {}, {}
        for e, f in itertools.combinations(range(E), 2):
            Ce, Cf = clades[e], clades[f]
            if not Ce & Cf: dis[Ce | Cf] = dis.get(Ce | Cf, 0.0) + pair_w(e, f)
            elif Cf < Ce or Ce < Cf:
                D = (Ce - Cf) if Cf < Ce else (Cf - Ce); nest[D] = nest.get(D, 0.0) + pair_w(e, f) * (1 - al * al)
        for s in all_splits(n):
            if s in T: continue
            C, D = s, X - s; taus += 1
            u = LU.get(s, 0.0); rc, rdd = LR.get(C, 0.0), LR.get(D, 0.0)
            ok &= rc <= u + 1e-15 and rdd <= u + 1e-15
            worst = min(worst, u - max(rc, rdd)); strict += u > max(rc, rdd) + 1e-15
            pu = dis.get(C, 0.0) + dis.get(D, 0.0) + nest.get(C, 0.0) + nest.get(D, 0.0)
            decomp_err = max(decomp_err, abs(pu - u), abs(al * al * dis.get(C, 0.0) + nest.get(C, 0.0) - rc),
                             abs(al * al * dis.get(D, 0.0) + nest.get(D, 0.0) - rdd))
    ok &= outcome_fail == 0 and decomp_err < 1e-14 and strict > 0
    return dict(passed=bool(ok), numbers=dict(trees=trees, false_splits=taus, strict=strict, smallest_slack=worst,
                                             outcome_failures=outcome_fail, decomposition_max_error=decomp_err),
                instances=f'{trees} random rooted trees ($n=5$--$7$, a third with polytomies), random rates, $\\alpha\\in\\{{0,1\\}}$ or '
                          f'random; all outcomes with at most two innovations; all {taus} splits not of the tree',
                result=f'second-order rooted support of each cluster $\\le$ unrooted support of the split in all (strict in {strict}); '
                       f'per outcome, a rooted proposal of $C$ or $D$ implies the unrooted proposal of $\\tau$; the decomposition of '
                       f'the proof (disjoint pairs $\\alpha^2$ against $1$, nested pairs $1-\\alpha^2$ in both) equals the enumeration '
                       f'(largest difference ${tex_e(decomp_err) if decomp_err else "0"}$)')


CHECKS = dict(exclusive=check_exclusive, nphard=check_nphard, buneman=check_buneman, threshold=check_threshold,
              pairwise=check_pairwise, consistency=check_consistency, samplesize=check_samplesize, blocks=check_blocks,
              gammablocks=check_gammablocks, textillegible=check_textillegible, bediereffect=check_bediereffect,
              bedier=check_bedier, rival=check_rival, dimfree=check_dimfree,
              dominance=check_dominance, modelglobal=check_modelglobal, lba=check_lba)


def run(name):
    t0 = time.time()
    try:
        r = CHECKS[name]()
    except Exception as e:                                   # a crash is a failure
        r = dict(passed=False, instances='', result=f'exception: {e!r}', numbers={})
    r['seconds'] = round(time.time() - t0, 1)
    print(name, 'PASS' if r['passed'] else 'FAIL', f"({r['seconds']} s):", r['result'], flush=True)
    return r


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', 'data', 'audit_new_checks.json')
    names = sys.argv[2:] or list(CHECKS)
    res = {k: run(k) for k in names}
    json.dump(res, open(out, 'w'), indent=1, default=float)
    bad = [k for k, v in res.items() if not v['passed']]
    print('ALL NEW CHECKS PASS' if not bad else f'FAILED: {bad}'); sys.exit(1 if bad else 0)
