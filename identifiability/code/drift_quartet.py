"""Root drift on the rooted quartet (((1,3),2),4), root at the vertex joining {1,2,3} and 4, attraction model.
Edges: pendants p1..p4, q above {1,3}, s above {1,2,3}. Exact supports of the true and competing splits/clusters by
summing over all outcomes (sympy), and their leading orders with q = rho eps^2, s = sigma eps, p3 = b eps,
p4 = tau eps, p1 = p2 = 0. Also a Monte Carlo check at a random point. Usage: drift_quartet.py OUT.json"""
import sys, json, itertools, random
import sympy as sp

al, q, s_, p1, p2, p3, p4, eps, rho, sig, b, tau = sp.symbols('alpha q s p1 p2 p3 p4 epsilon rho sigma b tau', nonnegative=True)
X = frozenset({1, 2, 3, 4})
edges = [(frozenset({1}), p1), (frozenset({2}), p2), (frozenset({3}), p3), (frozenset({4}), p4),
         (frozenset({1, 3}), q), (frozenset({1, 2, 3}), s_)]


def classes(out):
    val = {}
    for leaf in X:
        v = 'o'
        for k in sorted(range(len(edges)), key=lambda k: -len(edges[k][0])):
            if leaf in edges[k][0] and out[k] != 'none': v = 'a' if out[k] == 'a' else ('f', k)
        val[leaf] = v
    C = {}
    for l, v in val.items(): C.setdefault(v, set()).add(l)
    return C


def props(out):
    C = classes(out); U, R = set(), set()
    for v, c in C.items():
        c = frozenset(c)
        if len(c) == 2: U.add(c if 1 in c else X - c)
        if v != 'o' and 2 <= len(c) <= 3: R.add(c)
    return U, R


mU, mR = {}, {}
for out in itertools.product(['none', 'a', 'f'], repeat=6):
    pr = sp.Integer(1)
    for (c, rate), o in zip(edges, out): pr *= (1 - rate) if o == 'none' else rate * (al if o == 'a' else 1 - al)
    U, R = props(out)
    for x in U: mU[x] = mU.get(x, 0) + pr
    for x in R: mR[x] = mR.get(x, 0) + pr
sc = {q: rho * eps ** 2, s_: sig * eps, p3: b * eps, p4: tau * eps, p1: 0, p2: 0}
nm = lambda x: ''.join(map(str, sorted(x)))
res = {}
for tag, D in (('unrooted', mU), ('rooted', mR)):
    res[tag] = {}
    for x, e in sorted(D.items(), key=lambda kv: nm(kv[0])):
        lead = sp.factor(sp.series(sp.expand(e.subs(sc)), eps, 0, 3).removeO())
        res[tag][nm(x)] = str(lead); print(tag, nm(x), lead)
random.seed(3); vals = {al: 0.4, q: 0.05, s_: 0.2, p1: 0.03, p2: 0.06, p3: 0.35, p4: 0.1}; N = 300000; cU, cR = {}, {}
for _ in range(N):
    out = []
    for c, rate in edges:
        x = random.random(); out.append('none' if x >= float(vals[rate]) else ('a' if random.random() < float(vals[al]) else 'f'))
    U, R = props(out)
    for x in U: cU[x] = cU.get(x, 0) + 1
    for x in R: cR[x] = cR.get(x, 0) + 1
z = max(abs(float(e.subs(vals)) - C.get(x, 0) / N) / ((float(e.subs(vals)) * (1 - float(e.subs(vals))) / N) ** 0.5 + 1e-9)
        for D, C in ((mU, cU), (mR, cR)) for x, e in D.items())
res['monte_carlo_max_z'] = round(z, 2); print('MC max z', round(z, 2))
json.dump(res, open(sys.argv[1], 'w'), indent=1)
