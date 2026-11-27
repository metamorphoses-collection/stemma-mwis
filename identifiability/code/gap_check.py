"""Largest-gap threshold. Deterministic lemma checked on random instances: let true values lie in [b, 1], false in
[0, a] (0 included as an anchor), with b - a = G larger than every other gap g2 between consecutive population values;
if every empirical value is within eta < (G - g2)/4 of its population value, the largest gap of the empirical values
(with 0) lies between the largest false and the smallest true value, so a cut at its midpoint keeps exactly the true
ones. Also applies the rule to the Lai departures. Exits non-zero on failure. Usage: gap_check.py OUT.json"""
import sys, json, random, os
from collections import Counter
random.seed(9); ok = True; n_inst = 0
for _ in range(20000):
    k, f = random.randint(1, 8), random.randint(0, 12)
    tv = [random.uniform(0.2, 1) for _ in range(k)]; fv = [random.uniform(0, 0.2) for _ in range(f)]
    vals = sorted([0.0] + fv + tv)
    a, b = max([0.0] + fv), min(tv)
    gaps = [(y - x, x, y) for x, y in zip(vals, vals[1:])]
    G = b - a; g2 = max([g for g, x, y in gaps if not (x == a and y == b)], default=0.0)
    if G <= g2: continue
    eta = random.uniform(0, (G - g2) / 4 * 0.999)
    et = [v + random.uniform(-eta, eta) for v in tv]; ef = [max(0.0, v + random.uniform(-eta, eta)) for v in fv]
    ev = sorted([0.0] + ef + et)
    g, x, y = max((y - x, x, y) for x, y in zip(ev, ev[1:]))
    c = (x + y) / 2
    ok &= all(v > c for v in et) and all(v <= c for v in ef); n_inst += 1
# the Lai: departures counted against the printed text, all candidate clusters (unrooted sides)
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'lai_ombre', 'units_v2.json')
U = json.load(open(P)); W = 'ABCDEFG'; cnt = Counter()
for u in U:
    for kk, g in enumerate(u['groups']):
        g = ''.join(x for x in g if x in W)
        if kk != u['lemma_group'] and 2 <= len(g) <= 5:
            s = frozenset(g); s = s if 'A' not in s else frozenset(W) - s; cnt[s] += 1
L = len(U); sup = sorted(cnt.values(), reverse=True)
vals = sorted([0] + list(cnt.values()))
g, x, y = max((y - x, x, y) for x, y in zip(vals, vals[1:]))
kept = sorted((''.join(sorted(s)), v) for s, v in cnt.items() if v > (x + y) / 2)
res = dict(lemma_passed=bool(ok), instances=n_inst, lai_top_supports=sup[:8], lai_largest_gap=[x, y], lai_kept_by_gap_rule=kept)
print(res)
json.dump(res, open(sys.argv[1], 'w'), indent=1)
sys.exit(0 if ok else 1)
