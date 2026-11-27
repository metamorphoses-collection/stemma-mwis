"""Attraction model on the quartet: exact expected supports under the unrooted (reading-vs-rest) and rooted
(departures from the original reading) codings, and their leading orders in the long-branch scaling.
Model: each edge independently innovates with its rate; an innovation adopts the attractive reading 'a' with
probability alpha, otherwise a reading of its own; a witness reads the value of the lowest innovating edge above it,
or the original 'o' if none. alpha = 1 is the irreversible model of Section 7.
Exits non-zero unless both alpha = 1 reductions hold exactly and the Monte Carlo
largest |z| is below 5. Usage: kappa_quartet.py OUT.json"""
import sys, json, itertools
import sympy as sp

al, r1, r2, p1, p2, p3, p4, eps, a, b, rho = sp.symbols('alpha r1 r2 p1 p2 p3 p4 epsilon a b rho', nonnegative=True)
edges = [('p1', p1, (1,)), ('p2', p2, (2,)), ('p3', p3, (3,)), ('p4', p4, (4,)), ('h1', r1, (1, 2)), ('h2', r2, (3, 4))]


def leaf_values(out):
    val = {}
    for leaf in (1, 2, 3, 4):
        v = 'o'
        half = 'h1' if leaf in (1, 2) else 'h2'
        if out[half] != 'none': v = 'a' if out[half] == 'attract' else 'f_' + half
        if out[f'p{leaf}'] != 'none': v = 'a' if out[f'p{leaf}'] == 'attract' else f'f_p{leaf}'
        val[leaf] = v
    return val


def supports():
    U, R = {}, {}
    for outs in itertools.product(['none', 'attract', 'fresh'], repeat=6):
        out = dict(zip([e[0] for e in edges], outs))
        pr = sp.Integer(1)
        for (name, rate, _), o in zip(edges, outs):
            pr *= (1 - rate) if o == 'none' else rate * (al if o == 'attract' else 1 - al)
        val = leaf_values(out)
        classes = {}
        for leaf, v in val.items(): classes.setdefault(v, set()).add(leaf)
        uns, roo = set(), set()
        for v, C in classes.items():
            if len(C) == 2:
                uns.add(frozenset(C) if 1 in C else frozenset({1, 2, 3, 4}) - frozenset(C))
            if v != 'o' and 2 <= len(C) <= 3: roo.add(frozenset(C))
        for s in uns: U[s] = U.get(s, 0) + pr
        for c in roo: R[c] = R.get(c, 0) + pr
    return U, R


U, R = supports()
name = lambda s: ''.join(map(str, sorted(s)))
lb = {p2: 0, p4: 0, r1: rho * eps ** 2, r2: rho * eps ** 2, p1: a * eps, p3: b * eps}
res = dict(unrooted={}, rooted={})
for tag, D in (('unrooted', U), ('rooted', R)):
    for s, e in sorted(D.items(), key=lambda kv: name(kv[0])):
        full = sp.factor(sp.simplify(e))
        sub = sp.factor(sp.simplify(e.subs({p2: 0, p4: 0, r2: r1})))
        lead = sp.series(sp.expand(e.subs(lb)), eps, 0, 4).removeO()
        res[tag][name(s)] = dict(exact=str(full), long_branch_config=str(sub), leading=str(sp.factor(sp.expand(lead))))
        print(tag, name(s), '| p2=p4=0, r1=r2:', sub, '| leading:', sp.factor(sp.expand(lead)))
# alpha = 1 check against the irreversible model (Proposition quartet)
irr12 = (1 - r2) * (1 - p3) * (1 - p4) * (r1 + (1 - r1) * p1 * p2) + (1 - r1) * (1 - p1) * (1 - p2) * (r2 + (1 - r2) * p3 * p4)
irr13 = (1 - r1) * (1 - r2) * (p1 * p3 * (1 - p2) * (1 - p4) + p2 * p4 * (1 - p1) * (1 - p3))
chk = {'unrooted_12|34_alpha1': sp.simplify(U[frozenset({1, 2})].subs(al, 1) - irr12) == 0,
       'unrooted_13|24_alpha1': sp.simplify(U[frozenset({1, 3})].subs(al, 1) - irr13) == 0}
# Monte Carlo check of every exact support at a random parameter point (can fail)
import random
random.seed(7)
vals = {al: 0.3, r1: 0.07, r2: 0.11, p1: 0.35, p2: 0.05, p3: 0.4, p4: 0.08}
N = 400000; cU, cR = {}, {}
for _ in range(N):
    out = {}
    for name_, rate, _ in edges:
        x = random.random(); q = float(vals[rate])
        out[name_] = 'none' if x >= q else ('attract' if random.random() < float(vals[al]) else 'fresh')
    val = leaf_values(out); classes = {}
    for leaf, v in val.items(): classes.setdefault(v, set()).add(leaf)
    ks = set()
    for v, C in classes.items():
        if len(C) == 2: ks.add(frozenset(C) if 1 in C else frozenset({1, 2, 3, 4}) - frozenset(C))
        if v != 'o' and 2 <= len(C) <= 3: cR[frozenset(C)] = cR.get(frozenset(C), 0) + 1
    for k in ks: cU[k] = cU.get(k, 0) + 1          # a split is counted once per unit
worst = 0.0
for D, C in ((U, cU), (R, cR)):
    for k, e in D.items():
        ex = float(e.subs(vals)); mc = C.get(k, 0) / N; se = (ex * (1 - ex) / N) ** 0.5 + 1e-9
        worst = max(worst, abs(ex - mc) / se)
chk['monte_carlo_max_z'] = round(worst, 2); chk['monte_carlo_ok'] = worst < 5
print(chk)
res['checks'] = {k: (bool(v) if isinstance(v, bool) else v) for k, v in chk.items()}
passed = bool(chk['unrooted_12|34_alpha1'] and chk['unrooted_13|24_alpha1'] and chk['monte_carlo_ok'])
res['passed'] = passed
print('ALL PASS' if passed else 'FAIL', '(alpha = 1 reduction to the irreversible model; Monte Carlo |z| < 5)')
json.dump(res, open(sys.argv[1], 'w'), indent=1)
sys.exit(0 if passed else 1)
