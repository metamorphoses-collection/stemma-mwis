r"""Machine check of the three-mechanism theorem for the attraction model on random rooted binary trees (n = 5, 6).
For every outcome with at most two marked edges, the proposed splits (unrooted coding) and clusters (rooted coding)
are compared with the classification of the theorem:
  (i)   at most one mark: everything proposed is true or trivial;
  (ii)  two marks on edges with disjoint clades: unrooted proposes (C_e u C_f)|rest whatever the readings; rooted
        proposes the cluster C_e u C_f iff both readings are the attractive one;
  (iii) two marks with C_f strictly inside C_e: both codings propose C_e \ C_f unless both are attractive.
Then, for every split/cluster tau, exact m_tau (all 3^|E| outcomes) is checked to lie between the two-mark formula
L_tau = sum over pairs producing tau of pi_e pi_f w(e,f) prod_{g not e,f} (1-pi_g) (+ the one-mark term for true
ones) and L_tau + P(at least three marks). Exits non-zero on failure. Usage: mechanisms_check.py OUT.json"""
import sys, json, itertools, random, math

random.seed(11)


def rooted_tree(n):
    nodes = [frozenset({i}) for i in range(n)]; cl = []
    while len(nodes) > 2:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; cl.append(c); nodes.append(c)
    return [frozenset({i}) for i in range(n)] + cl + nodes   # root = vertex joining the last two; its two child edges included


def propose(n, clades, outcome):
    X = frozenset(range(n)); val = {}
    order = sorted(range(len(clades)), key=lambda k: -len(clades[k]))      # higher edges first
    for leaf in range(n):
        v = 'o'
        for k in order:
            if leaf in clades[k] and outcome[k] != 'none':
                v = 'a' if outcome[k] == 'a' else ('f', k)                     # lowest marked edge wins (later = lower)
        val[leaf] = v
    classes = {}
    for leaf, v in val.items(): classes.setdefault(v, set()).add(leaf)
    U, Rt = set(), set()
    for v, C in classes.items():
        C = frozenset(C)
        if 2 <= len(C) <= n - 2: U.add(C if 0 not in C else X - C)
        if v != 'o' and 2 <= len(C) <= n - 1: Rt.add(C)
    return U, Rt


ok = True; report = []
for trial in range(40):
    n = random.choice([5, 6]); clades = rooted_tree(n); X = frozenset(range(n)); E = len(clades)
    canon = lambda C: C if 0 not in C else X - C
    trueU = {canon(c) for c in clades if 2 <= len(c) <= n - 2}
    trueR = {c for c in clades if 2 <= len(c) <= n - 1}
    is_clade = lambda C: C in set(clades)
    rates = [random.uniform(0.02, 0.3) for _ in clades]; al = random.uniform(0, 1)
    # (i)-(iii): deterministic classification on <=2 marks
    for k in range(E):
        for v in ('a', 'f'):
            out = ['none'] * E; out[k] = v; U, Rt = propose(n, clades, out)
            ok &= U <= trueU and Rt <= trueR
    for e, f in itertools.combinations(range(E), 2):
        Ce, Cf = clades[e], clades[f]
        for ve, vf in itertools.product('af', repeat=2):
            out = ['none'] * E; out[e] = ve; out[f] = vf; U, Rt = propose(n, clades, out)
            if not (Ce & Cf):
                Un = Ce | Cf
                if ve == vf == 'a':                       # one class C_e u C_f, and the original class
                    expU = {canon(Un)} if 2 <= len(Un) <= n - 2 else set()
                else:
                    expU = {canon(Ce), canon(Cf)} | ({canon(Un)} if 2 <= len(Un) <= n - 2 else set())
                    expU = {s for s in expU if 2 <= len(s) <= n - 2}
                expR = {c for c in (Ce, Cf) if 2 <= len(c) <= n - 1} if not (ve == vf == 'a') else ({Un} if 2 <= len(Un) <= n - 1 else set())
                if not (U == expU and Rt == expR): print('disjoint mismatch', sorted(map(sorted,U)), sorted(map(sorted,expU)), sorted(map(sorted,Rt)), sorted(map(sorted,expR)), ve, vf); ok = False
            else:
                big, small = (Ce, Cf) if Cf < Ce else (Cf, Ce)
                D = big - small
                if ve == vf == 'a':
                    expU = {canon(big)} if 2 <= len(big) <= n - 2 else set(); expR = {big} if 2 <= len(big) <= n - 1 else set()
                else:
                    expU = {canon(s) for s in (D, small) if 2 <= len(s) <= n - 2} | ({canon(big)} if 2 <= len(big) <= n - 2 else set())
                    expR = {s for s in (D, small) if 2 <= len(s) <= n - 1}
                ok &= U == expU and Rt == expR
    # exact supports vs two-mark bounds
    mU, mR = {}, {}
    P3 = 0.0
    for out in itertools.product(['none', 'a', 'f'], repeat=E):
        pr = 1.0
        for k, o in enumerate(out): pr *= (1 - rates[k]) if o == 'none' else rates[k] * (al if o == 'a' else 1 - al)
        if sum(o != 'none' for o in out) >= 3: P3 += pr
        U, Rt = propose(n, clades, list(out))
        for s in U: mU[s] = mU.get(s, 0) + pr
        for c in Rt: mR[c] = mR.get(c, 0) + pr
    LU, LR = {}, {}
    for e, f in itertools.combinations(range(E), 2):
        base = rates[e] * rates[f] * math.prod(1 - rates[g] for g in range(E) if g not in (e, f))
        for ve, vf in itertools.product('af', repeat=2):
            w = (al if ve == 'a' else 1 - al) * (al if vf == 'a' else 1 - al)
            out = ['none'] * E; out[e] = ve; out[f] = vf; U, Rt = propose(n, clades, out)
            for s in U: LU[s] = LU.get(s, 0) + base * w
            for c in Rt: LR[c] = LR.get(c, 0) + base * w
    for k in range(E):
        base = rates[k] * math.prod(1 - rates[g] for g in range(E) if g != k)
        for v in 'af':
            w = al if v == 'a' else 1 - al
            out = ['none'] * E; out[k] = v; U, Rt = propose(n, clades, out)
            for s in U: LU[s] = LU.get(s, 0) + base * w
            for c in Rt: LR[c] = LR.get(c, 0) + base * w
    for M, L in ((mU, LU), (mR, LR)):
        for s in set(M) | set(L):
            ok &= L.get(s, 0) - 1e-12 <= M.get(s, 0) <= L.get(s, 0) + P3 + 1e-12
    report.append(dict(n=n, edges=E, alpha=round(al, 3), P3=P3))
print('ALL PASS' if ok else 'FAIL', len(report), 'trees')
json.dump(dict(passed=bool(ok), trees=report), open(sys.argv[1], 'w'), indent=1)
sys.exit(0 if ok else 1)
