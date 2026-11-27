r"""Check of a claimed dimension-free bound (Theo, 22 Nov 2026) for the cutting-plane estimator under independent
illegibility: for a binary stemma T' and a non-trivial full split sigma' = A|B not in Sigma(T'), with each witness
illegible independently with probability delta,
   P( sigma'|O is non-trivial and displayed by T' ) <= 2 delta   for delta <= 1/2.
Exact enumeration over all legibility patterns, all false splits, random binary trees, n = 5..8.
Usage: tolerance_check.py OUT.json"""
import sys, json, itertools, random
import networkx as nx
random.seed(20261122)
def random_binary_tree(n):
    T = nx.Graph(); T.add_edges_from([(0, 'r'), (1, 'r'), (2, 'r')]); nxt = 0
    for leaf in range(3, n):
        u, v = random.choice(list(T.edges())); m = ('i', nxt); nxt += 1
        T.remove_edge(u, v); T.add_edges_from([(u, m), (m, v), (m, leaf)])
    return T
def splits(T, n):
    S = set()
    for u, v in T.edges():
        H = T.copy(); H.remove_edge(u, v)
        side = frozenset(x for x in nx.node_connected_component(H, u) if isinstance(x, int))
        if 2 <= len(side) <= n - 2: S.add(min(side, frozenset(range(n)) - side, key=sorted))
    return S
def displayed(T, splitsT, A, B):
    # A|B (partial) displayed iff some split of T restricts to it
    for C in splitsT:
        D = frozenset(range(n)) - C
        if (A <= C and B <= D) or (A <= D and B <= C): return True
    return False
out = {}; worst = (0, None)
for n in (5, 6, 7, 8):
    X = frozenset(range(n)); res = []
    for trial in range(25 if n < 8 else 10):
        T = random_binary_tree(n); ST = splits(T, n)
        STfull = ST | {X - s for s in ST}
        allsides = [frozenset(c) for k in range(2, n - 1) for c in itertools.combinations(range(n), k) if 0 in c]
        for A in allsides:
            B = X - A
            if A in ST or B in ST: continue
            # probability as polynomial evaluated at several deltas
            pats = []
            for r in range(n + 1):
                for ill in itertools.combinations(range(n), r):
                    O = X - frozenset(ill); a, b = A & O, B & O
                    if len(a) >= 2 and len(b) >= 2 and displayed(T, STfull, a, b): pats.append(r)
            for d in (0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
                p = sum(d ** r * (1 - d) ** (n - r) for r in pats)
                ratio = p / (2 * d)
                if ratio > worst[0]: worst = (ratio, dict(n=n, A=sorted(A), delta=d, p=p))
    out[n] = 'done'
print('worst P/(2 delta):', worst)
json.dump(dict(worst_ratio=worst[0], worst_case=worst[1]), open(sys.argv[1], 'w'), indent=1)
