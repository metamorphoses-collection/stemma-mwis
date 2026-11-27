r"""The trio D, E, F of the Lai under an unequal-rate null (prop:polytomynull with A in place of alpha^2).
At the vertex joining D, E, F (not the root), m(DE) = A q_D q_E + (1-A) q_DEF q_F, and symmetrically, to O(eps^3).
Each q is read from the number of places at which exactly that set departs from the printed text. For a grid of A,
the exact multinomial test of the observed counts (DE, EF, DF) against the null probabilities is computed (sum of
the probabilities of all outcomes no more likely than the observed one). Also: the equal-rate selection-corrected
p-value after merging adjacent units and after halving the counts. Usage: rival_unequal.py OUT.json"""
import sys, json, os, math
from math import lgamma, exp, log
from collections import Counter
W = 'ABCDEFG'; X = frozenset(W)
U = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'lai_ombre', 'units_v2.json')))
exact = Counter(); pos = {'DE': [], 'EF': [], 'DF': []}
for i, u in enumerate(U):
    for k, g in enumerate(u['groups']):
        g = ''.join(x for x in g if x in W)
        if k == u['lemma_group'] or not g: continue
        exact[''.join(sorted(g))] += 1
        if 2 <= len(g) <= 5:
            s = frozenset(g); s = s if 'A' not in s else X - s; key = ''.join(sorted(s))
            if key in pos: pos[key].append(i)
obs = tuple(len(pos[k]) for k in ('DE', 'EF', 'DF')); N = sum(obs)
q = {k: exact[k] for k in ('D', 'E', 'F', 'DEF')}


def lp(x, pr): return lgamma(N + 1) - sum(lgamma(v + 1) for v in x) + sum(v * log(p) for v, p in zip(x, pr))


def gof(pr):
    lo = lp(obs, pr); return sum(exp(l) for a in range(N + 1) for b in range(N + 1 - a)
                                 for l in [lp((a, b, N - a - b), pr)] if l <= lo + 1e-12)


def tail(a, n): return sum(math.comb(n, i) for i in range(a, n + 1)) / 2 ** n


def tri(x):                      # equal-rate selection-corrected statistic, as in the trinomial test of the paper
    M = sum(x); s = sorted(x, reverse=True); so = tail(s[0], s[0] + s[1]); p = 0
    for a in range(M + 1):
        for b in range(M + 1 - a):
            z = sorted((a, b, M - a - b), reverse=True)
            if tail(z[0], z[0] + z[1]) <= so + 1e-15:
                p += math.factorial(M) / (math.factorial(a) * math.factorial(b) * math.factorial(M - a - b)) / 3 ** M
    return p


rows = []
for A in (0, 0.02, 0.05, 0.1, 0.2, 0.34, 0.39, 0.5, 1.0):
    m = (A * q['D'] * q['E'] + (1 - A) * q['DEF'] * q['F'], A * q['E'] * q['F'] + (1 - A) * q['DEF'] * q['D'],
         A * q['D'] * q['F'] + (1 - A) * q['DEF'] * q['E'])
    pr = [v / sum(m) for v in m]
    rows.append(dict(A=A, expected=[round(N * v, 2) for v in pr], heaviest=['DE', 'EF', 'DF'][pr.index(max(pr))], p=gof(pr)))
adj = {k: sum(1 for a, b in zip(v, v[1:]) if b - a <= 1) for k, v in pos.items()}
merged = tuple(len(pos[k]) - adj[k] for k in ('DE', 'EF', 'DF'))
lines = len({(U[i]['page'], U[i]['line']) for i in pos['DE']})
res = dict(observed=dict(zip(('DE', 'EF', 'DF'), obs)), exact_departures=q, null=rows,
           equal_rate_p=tri(obs), DE_distinct_lines=lines, adjacent_pairs=adj, merged=merged, merged_p=tri(merged),
           halved_p=[tri((20, 11, 6)), tri((21, 12, 7))])
json.dump(res, open(sys.argv[1], 'w'), indent=1); print(json.dumps(res, indent=1))
