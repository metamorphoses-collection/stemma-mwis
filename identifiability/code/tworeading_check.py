r"""Machine check of the two-reading coding (prop:tworeading).
A unit proposes a split A|B only when A and B are exactly its two reading classes.
(i) On random rooted binary trees (n = 5, 6, 7), every outcome with at most two innovations (none / attractive /
    own reading on each innovating edge) is enumerated: one innovation proposes a true split or nothing; two
    innovations propose a split of no edge only for a disjoint pair with both readings attractive, and that split
    is (C_e u C_f)|rest; a nested pair proposes nothing or the split of e (both attractive).
(ii) On the quartet 12|34 (p2 = p4 = 0, halves r), the exact supports summed over the 81 edge outcomes equal
    m(13|24) = alpha^2 p1 p3 (1-r)^2, m(14|23) = 0, and m(12|34) = 2r + O(eps^3) under r = rho eps^2, p = a eps, b eps.
Exits non-zero on failure. Usage: tworeading_check.py OUT.json"""
import sys, json, itertools, random
import sympy as sp
random.seed(28)


def rooted_tree(n):
    nodes = [frozenset({i}) for i in range(n)]; cl = []
    while len(nodes) > 2:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; cl.append(c); nodes.append(c)
    return [frozenset({i}) for i in range(n)] + cl          # every edge of the rooted tree: singletons and clades below the root


def classes(n, clades, out):
    order = sorted(range(len(clades)), key=lambda k: -len(clades[k])); cls = {}
    for leaf in range(n):
        v = 'o'
        for k in order:
            if leaf in clades[k] and out.get(k) is not None:
                v = 'A' if out[k] == 'a' else ('f', k)
        cls.setdefault(v, set()).add(leaf)
    return [frozenset(c) for c in cls.values()]


def canon(A, n):
    X = frozenset(range(n)); return min(A, X - A, key=lambda z: sorted(z))


def two_reading(n, cls):
    if len(cls) == 2 and all(len(c) >= 2 for c in cls): return {canon(cls[0], n)}
    return set()


ok = True; trials = []; coincident_false = 0
for n in (5, 6, 7):
    for _ in range(8):
        clades = rooted_tree(n)          # edges: every non-root clade (singletons and internal)
        # edges of the rooted tree: all clades except the two root children merge into X; include the two root children
        X = frozenset(range(n)); edges = [c for c in clades if c != X]
        true = {canon(c, n) for c in edges if 2 <= len(c) <= n - 2}
        bad = 0
        for k in range(len(edges)):
            for o in ('a', 'f'):
                P = two_reading(n, classes(n, edges, {k: o}))
                bad += any(s not in true for s in P)
        for e, f in itertools.combinations(range(len(edges)), 2):
            Ce, Cf = edges[e], edges[f]
            for oe, of in itertools.product('af', repeat=2):
                P = two_reading(n, classes(n, edges, {e: oe, f: of}))
                false = [s for s in P if s not in true]
                if not (Ce & Cf):
                    want = set()
                    if oe == of == 'a' and len(Ce | Cf) <= n - 2 and len(Ce | Cf) >= 2:
                        u = canon(Ce | Cf, n); want = {u} if u not in true else set(); coincident_false += bool(want)
                    bad += set(false) != want
                else:
                    bad += bool(false)
        trials.append(dict(n=n, edges=len(edges), failures=bad)); ok &= bad == 0

# (ii) quartet
r, p1, p3, al, rho, a, b, eps = sp.symbols('r p1 p3 alpha rho a b epsilon', positive=True)
rate = {'h12': r, 'h34': r, 'p1': p1, 'p2': 0, 'p3': p3, 'p4': 0}
above = {1: ['p1', 'h12'], 2: ['p2', 'h12'], 3: ['p3', 'h34'], 4: ['p4', 'h34']}
E = list(rate); m = {'12': 0, '13': 0, '14': 0}
for out in itertools.product([0, 1, 2], repeat=6):
    pr = 1
    for e, o in zip(E, out):
        q = rate[e]; pr *= (1 - q) if o == 0 else (q * al if o == 1 else q * (1 - al))
    if pr == 0: continue
    o = dict(zip(E, out)); rd = {}
    for w in (1, 2, 3, 4):
        v = 0
        for e in above[w]:
            if o[e] == 1: v = 'A'; break
            if o[e] == 2: v = e; break
        rd[w] = v
    if rd[1] == rd[2] and rd[3] == rd[4] and rd[1] != rd[3]: m['12'] += pr
    if rd[1] == rd[3] and rd[2] == rd[4] and rd[1] != rd[2]: m['13'] += pr
    if rd[1] == rd[4] and rd[2] == rd[3] and rd[1] != rd[2]: m['14'] += pr
q1 = sp.simplify(m['13'] - al**2 * p1 * p3 * (1 - r)**2) == 0
q2 = sp.simplify(m['14']) == 0
lead = sp.series(sp.expand(m['12'].subs({r: rho * eps**2, p1: a * eps, p3: b * eps})), eps, 0, 3).removeO()
q3 = sp.simplify(lead - 2 * rho * eps**2) == 0
ok &= q1 and q2 and q3
ok &= coincident_false > 0       # the false case must occur, or the check could not fail
res = dict(trials=trials, coincident_false_cases=coincident_false, quartet=dict(m13=q1, m14_zero=q2, m12_leading=q3), ok=bool(ok))
json.dump(res, open(sys.argv[1], 'w'), indent=1); print(res)
sys.exit(0 if ok else 1)
