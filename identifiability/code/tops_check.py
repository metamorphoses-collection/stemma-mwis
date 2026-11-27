r"""Machine check of the corollary on polytomies that survive (cor:tops).
For random non-negative cluster weights on n = 5 witnesses (many weights zero), every maximum-weight laminar set of
proper clusters (sizes 2..n-1) is found by exhaustive enumeration. For each optimum S and each vertex V (a member of S,
or the full set) with k >= 3 children (maximal members of S and singletons strictly inside V), no union of two
children has positive weight. Exits non-zero on failure. Usage: tops_check.py OUT.json"""
import sys, json, itertools, random
random.seed(29)
n = 5; X = frozenset(range(n))
CL = [frozenset(c) for k in range(2, n) for c in itertools.combinations(range(n), k)]
lam = lambda a, b: not (a & b) or a <= b or b <= a
FAM = [()]
for k in range(1, n - 1):
    FAM += [F for F in itertools.combinations(CL, k) if all(lam(a, b) for a, b in itertools.combinations(F, 2))]


def children(V, S):
    inside = [c for c in S if c < V] + [frozenset({i}) for i in V]
    return [c for c in inside if not any(c < d for d in inside)]


ok = True; checked = polytomies = 0
for t in range(400):
    w = {c: (random.randint(1, 9) if random.random() < 0.35 else 0) for c in CL}
    vals = [sum(w[c] for c in F) for F in FAM]; best = max(vals)
    for F, v in zip(FAM, vals):
        if v != best or any(w[c] == 0 for c in F): continue       # optima using only candidates
        checked += 1
        for V in list(F) + [X]:
            ch = children(V, F)
            if len(ch) >= 3:
                polytomies += 1
                if any(w.get(a | b, 0) > 0 for a, b in itertools.combinations(ch, 2)): ok = False
ok &= polytomies > 0
res = dict(n=n, weightings=400, optima=checked, polytomies=polytomies, ok=bool(ok))
json.dump(res, open(sys.argv[1], 'w'), indent=1); print(res); sys.exit(0 if ok else 1)
