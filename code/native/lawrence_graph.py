"""Lawrence prologue (BHL 4759) conflict graph for a native Rydberg encoding.
Families = sets of witnesses sharing a non-lemma reading (side without T, the complete base text), counted
once per unit; keep support >= MIN. Conflict = the two families cross (intersect, neither contains the other).
Outputs the graph and brute-force MIS / MWIS. Usage: lawrence_graph.py <matrix.num> <out.json> [MIN]"""
import sys, json, itertools as it
from collections import Counter
sys.path.insert(0, __file__.rsplit('/', 2)[0] + '/benchmark')
from stemma_mwis import load_num
wit, rows = load_num(sys.argv[1]); MIN = int(sys.argv[3]) if len(sys.argv) > 3 else 2
t = wit.index('T'); cnt = Counter()
for r in rows:
    obs = [i for i in range(len(wit)) if r[i] != '?']
    fams = set()
    for v in set(r[i] for i in obs):
        A = frozenset(wit[i] for i in obs if r[i] == v)
        if 'T' in A: A = frozenset(wit[i] for i in obs) - A
        if 2 <= len(A) <= len(obs) - 2 or (len(A) >= 2 and len(obs) - len(A) >= 1 and 'T' not in A and len(A) < len(wit) - 1):
            fams.add(A)
    for A in fams: cnt[A] += 1
F = sorted([A for A, c in cnt.items() if c >= MIN], key=lambda A: (-cnt[A], sorted(A)))
names = ['+'.join(sorted(A)) for A in F]; W = {names[i]: cnt[F[i]] for i in range(len(F))}
cross = lambda a, b: bool(a & b) and not (a <= b or b <= a)
edges = [(i, j) for i, j in it.combinations(range(len(F)), 2) if cross(F[i], F[j])]
n = len(F); adj = [[0] * n for _ in range(n)]
for i, j in edges: adj[i][j] = adj[j][i] = 1
ind = [s for k in range(n + 1) for s in it.combinations(range(n), k) if all(not adj[a][b] for a, b in it.combinations(s, 2))]
m = max(map(len, ind)); mis = [s for s in ind if len(s) == m]
best = max(ind, key=lambda s: sum(W[names[i]] for i in s))
print('families (support):', [(names[i], W[names[i]]) for i in range(n)])
print('conflicts:', [(names[i], names[j]) for i, j in edges])
print(f'n={n} edges={len(edges)}  MIS size {m}, {len(mis)} maximum sets:', [[names[i] for i in s] for s in mis])
print('MWIS:', [names[i] for i in best], 'weight', sum(W[names[i]] for i in best))
json.dump(dict(names=names, W=W, n=n, edges=edges, adj=adj), open(sys.argv[2], 'w'), indent=1)
