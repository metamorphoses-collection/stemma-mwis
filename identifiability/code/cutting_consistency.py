r"""Theorem cuttingconsistency, checked exactly on small instances.

(a) A family of partial splits is jointly completable iff some binary tree displays every member (p = (O, S) is
    displayed by T' if some split of T' restricts to S | O\S on O). Hence the maximum-weight jointly completable family
    has the weight max over binary trees T' of W(T') = sum of w(p) over observed p displayed by T'.
    Check: random small collations on n = 5, 6 with random illegibility; brute force over all subfamilies, joint
    completability by backtracking over completions, against the maximum over all binary trees.
(b, c) With M(T') = E[number of distinct partial proposals of a unit displayed by T'], m the complete-data supports,
    mu = E N and mu2 = E N(N-1)/2 (N = number of splits a complete unit proposes), for every binary T' != T:
      M(T) - M(T') >= m(T) - m(T') - delta (4 mu + 2 (n-3) mu + 2 mu2).
    Check: exact enumeration of all innovation states (attraction model, unrooted coding) and all legibility patterns
    on random rooted binary trees with n = 5, random rates, alpha and delta. Usage: cutting_consistency.py OUT.json"""
import sys, json, random, itertools
from collections import defaultdict
random.seed(29)


def nontrivial_splits(n):
    out = []
    for r in range(2, n - 1):
        for c in itertools.combinations(range(1, n), r):
            if r <= n - 2: out.append(frozenset(c))          # side without witness 0
    return out


def compatible(A, B, X):
    return not (A & B and A - B and B - A and X - (A | B))


def binary_trees(n):
    X = frozenset(range(n)); S = nontrivial_splits(n); out = []
    for comb in itertools.combinations(S, n - 3):
        if all(compatible(a, b, X) for a, b in itertools.combinations(comb, 2)): out.append(frozenset(comb))
    return X, out


def canon(O, S):
    """partial split on O with side S; canonical side = the one containing min(O)."""
    m = min(O); side = S if m in S else O - S
    return (O, frozenset(side))


def displays(tree, p, X):
    O, S = p
    for A in tree:
        a = A & O
        if a == S or (O - a) == S: return True
    return False


def completable(fam, X):
    """backtracking: choose for each partial split a full extension, pairwise compatible."""
    fam = list(fam); n = len(X)
    def exts(p):
        O, S = p; U = sorted(X - O); out = []
        for bits in itertools.product((0, 1), repeat=len(U)):
            A = set(S) | {u for u, b in zip(U, bits) if b}
            A = frozenset(A) if 0 not in A else X - frozenset(A)
            if 2 <= len(A) <= n - 2: out.append(A)
        return out
    E = [exts(p) for p in fam]
    def bt(k, chosen):
        if k == len(fam): return True
        for A in E[k]:
            if all(compatible(A, B, X) for B in chosen) and bt(k + 1, chosen + [A]): return True
        return False
    return bt(0, [])


# ---------- (a) equivalence ----------
eq_ok = True; eq_tests = 0
for n in (5, 6):
    X, trees = binary_trees(n)
    for _ in range(25 if n == 5 else 8):
        w = defaultdict(int)
        for u in range(random.randint(6, 10)):
            O = frozenset(x for x in range(n) if random.random() > 0.25)
            if len(O) < 4: continue
            S = frozenset(random.sample(sorted(O), random.randint(2, len(O) - 2)))
            w[canon(O, S)] += random.randint(1, 3)
        P = list(w)
        if len(P) > 11: P = P[:11]
        best_fam = 0
        for r in range(len(P), -1, -1):
            for sub in itertools.combinations(P, r):
                wt = sum(w[p] for p in sub)
                if wt > best_fam and completable(sub, X): best_fam = wt
        best_tree = max(sum(w[p] for p in P if displays(t, p, X)) for t in trees)
        eq_ok &= best_fam == best_tree; eq_tests += 1
print('(a) equivalence:', eq_ok, eq_tests, flush=True)


# ---------- (b, c) the bound ----------
def rooted_binary(n):
    nodes = [frozenset([i]) for i in range(n)]; parent = {}; clade = {}
    k = n
    while len(nodes) > 1:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; nodes.append(c)
    # rebuild as list of clades (edges) excluding root
    return None


def random_rooted_clades(n):
    nodes = [frozenset([i]) for i in range(n)]; cl = []
    while len(nodes) > 1:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; nodes.append(c); cl += [a, b]
    return cl                                   # every non-root clade = one edge


def exact(n, cl, pi, alpha, delta, X, trees):
    m = defaultdict(float); M = defaultdict(float); mu = mu2 = 0.0
    E = list(range(len(cl)))
    legs = []
    for pat in itertools.product((0, 1), repeat=n):
        O = frozenset(i for i in range(n) if pat[i]); pr = 1.0
        for i in range(n): pr *= (1 - delta) if pat[i] else delta
        legs.append((O, pr))
    for st in itertools.product((0, 1, 2), repeat=len(E)):
        pr = 1.0
        for e, s in zip(E, st): pr *= (1 - pi[e]) if s == 0 else pi[e] * (alpha if s == 1 else 1 - alpha)
        read = {}
        for x in range(n):
            best, r = None, 0
            for e in E:
                if st[e] and x in cl[e] and (best is None or len(cl[e]) < len(cl[best])): best = e
            read[x] = 0 if best is None else ('a' if st[best] == 1 else ('o', best))
        cls = defaultdict(set)
        for x, r in read.items(): cls[r].add(x)
        full = set()
        for C in cls.values():
            C = frozenset(C); A = C if 0 not in C else X - C
            if 2 <= len(A) <= n - 2: full.add(A)
        for A in full: m[A] += pr
        N = len(full); mu += pr * N; mu2 += pr * N * (N - 1) / 2
        for O, pl in legs:
            if len(O) < 4: continue
            props = set()
            for C in cls.values():
                S = frozenset(C) & O
                if 2 <= len(S) <= len(O) - 2: props.add(canon(O, S))
            for t_i, t in enumerate(trees):
                M[t_i] += pr * pl * sum(1 for p in props if displays(t, p, X))
    return m, M, mu, mu2


n = 5; X, trees = binary_trees(n); rows = []; bound_ok = True
for trial in range(30):
    cl = random_rooted_clades(n)
    pi = [random.uniform(0.01, 0.25) for _ in cl]; alpha = random.uniform(0, 1); delta = random.uniform(0.01, 0.3)
    m, M, mu, mu2 = exact(n, cl, pi, alpha, delta, X, trees)
    T = frozenset(A if 0 not in A else X - A for A in cl if 2 <= len(A) <= n - 2)
    tT = trees.index(T)
    mT = sum(m[A] for A in T); slack_min = float('inf'); argmaxM = max(range(len(trees)), key=lambda i: M[i])
    for i, t in enumerate(trees):
        if i == tT: continue
        lhs = M[tT] - M[i]; rhs = mT - sum(m[A] for A in t) - delta * (4 * mu + 2 * (n - 3) * mu + 2 * mu2)
        slack_min = min(slack_min, lhs - rhs)
    ok = slack_min >= -1e-12; bound_ok &= ok
    rows.append(dict(alpha=round(alpha, 3), delta=round(delta, 3), slack_min=slack_min, true_maximises_M=argmaxM == tT,
                     true_maximises_m=max(trees, key=lambda t: sum(m[A] for A in t)) == T, passed=ok))
print('(c) bound holds:', bound_ok, 'min slack', min(r['slack_min'] for r in rows), flush=True)
# ---------- regression (GPT review, finding 11): unique maximiser whose split is never observed in full ----------
# leaves 1..5 = indices 0..4; T has cherries 12 and 45; a complete unit proposes 45 (p = 2/3) or 34 (p = 1/3); only
# leaf 5 can be illegible, with probability r. Then T is the unique maximum of M although m(12|345) = 0.
X5, trees5 = binary_trees(5); reg_ok = True
T5 = frozenset({X5 - frozenset({0, 1}), frozenset({3, 4})})
for r in (0.1, 0.3, 0.6):
    M = [0.0] * len(trees5)
    for full, pf in ((frozenset({3, 4}), 2 / 3), (frozenset({2, 3}), 1 / 3)):
        for leg5, pl in ((True, 1 - r), (False, r)):
            O = X5 if leg5 else X5 - {4}
            S = full & O
            if 2 <= len(S) <= len(O) - 2:
                pp = canon(O, S)
                for i, t in enumerate(trees5): M[i] += pf * pl * displays(t, pp, X5)
    iT = trees5.index(T5); best = max(M)
    unique = M[iT] == best and sum(abs(x - best) < 1e-12 for x in M) == 1
    reg_ok &= unique
print('(b) regression: T unique maximum with an unobserved split:', reg_ok, flush=True)
out = dict(regression_finding11=bool(reg_ok), equivalence_pass=bool(eq_ok), equivalence_tests=eq_tests, bound_pass=bool(bound_ok), bound_instances=len(rows),
           true_max_M=sum(r['true_maximises_M'] for r in rows), true_max_m=sum(r['true_maximises_m'] for r in rows), rows=rows)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
print({k: out[k] for k in out if k != 'rows'})
sys.exit(0 if eq_ok and bound_ok and reg_ok else 1)
