"""Rival test: for each family of a selected stemma, its strongest rival in conflict with it and the exact one-sided
binomial p-value of w(sigma) against w(rival) given their sum (valid because conflicting families are never proposed
by the same unit, so given N = w(sigma) + w(rival) the count w(sigma) is Binomial(N, m_sigma/(m_sigma+m_rival)); the
null is m_sigma = m_rival). Applied to the Lai (departures from the printed text, Bedier 1913). Usage: rival_test.py OUT"""
import sys, json, os, math
from collections import Counter


def binom_tail(k, n):            # P(Bin(n, 1/2) >= k)
    return sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n


W = 'ABCDEFG'; X = frozenset(W)
U = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'lai_ombre', 'units_v2.json')))
cnt = Counter()
for u in U:
    for kk, g in enumerate(u['groups']):
        g = ''.join(x for x in g if x in W)
        if kk != u['lemma_group'] and 2 <= len(g) <= 5:
            s = frozenset(g); cnt[s if 'A' not in s else X - s] += 1
def conflict(A, B):
    Ac, Bc = X - A, X - B
    return bool(A & B and A & Bc and Ac & B and Ac & Bc)
tree = [frozenset(x) if 'A' not in x else X - frozenset(x) for x in ('AB', 'CG', 'DE', 'DEF')]
out = []
for s in tree:
    rivals = sorted(((v, r) for r, v in cnt.items() if conflict(s, r)), key=lambda z: -z[0])
    v, r = rivals[0]; ws = cnt[s]
    out.append(dict(family=''.join(sorted(s)), support=ws, strongest_rival=''.join(sorted(r)), rival_support=v,
                    p_one_sided=binom_tail(ws, ws + v), crossing_rivals=len(rivals)))
    print(out[-1])
json.dump(out, open(sys.argv[1], 'w'), indent=1)
