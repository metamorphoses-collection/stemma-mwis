"""Bootstrap distribution over histories on a fixed candidate-family graph.
Resample the units with replacement B times; recompute each family's support (units whose split gives that family,
as in lawrence_graph.py); take the maximum-weight independent set among the fixed candidates (ties split equally).
Usage: bootstrap_histories.py <matrix.num> <graph.json> [B] -> prints P(history)"""
import sys, json, random, itertools as it
from collections import Counter
sys.path.insert(0, __file__.rsplit('/', 2)[0] + '/benchmark')
from stemma_mwis import load_num
wit, rows = load_num(sys.argv[1]); G = json.load(open(sys.argv[2])); B = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
names = G['names']; n = len(names); adj = G['adj']; fam = [frozenset(x.split('+')) for x in names]
def unit_fams(r):
    obs = [i for i in range(len(wit)) if r[i] != '?']; out = set()
    for v in set(r[i] for i in obs):
        A = frozenset(wit[i] for i in obs if r[i] == v)
        if 'T' in A: A = frozenset(wit[i] for i in obs) - A
        if A in fam: out.add(fam.index(A))
    return out
U = [unit_fams(r) for r in rows]
ind = [s for k in range(n + 1) for s in it.combinations(range(n), k) if all(not adj[a][b] for a, b in it.combinations(s, 2))]
maximal = [s for s in ind if not any(set(s) < set(t) for t in ind)]
lab = lambda s: '{' + ', '.join(names[i] for i in s) + '}'
full = Counter(i for u in U for i in u)
print('weights on full data:', {names[i]: full[i] for i in range(n)})
rng = random.Random(1); P = Counter()
for _ in range(B):
    w = Counter(i for u in (U[rng.randrange(len(U))] for _ in U) for i in u)
    val = {s: sum(w[i] for i in s) for s in maximal}; m = max(val.values()); top = [s for s in maximal if val[s] == m]
    for s in top: P[s] += 1 / len(top)
res = {lab(s): round(P[s] / B, 4) for s in maximal}
print(json.dumps(res, indent=1))
json.dump(dict(B=B, maximal=[lab(s) for s in maximal], P=res), open(sys.argv[2].replace('.json', '_bootstrap.json'), 'w'), indent=1)
