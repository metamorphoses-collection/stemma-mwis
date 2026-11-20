r"""Machine check of the fragmentary-witness theorem. Complete-data units follow the attraction model on a random
rooted binary tree (n = 5); each witness is then illegible independently with probability delta_i <= delta <= 1/2,
independently of the readings. The restriction-weighted support of a full split sigma is the probability that some
reading class, restricted to the legible witnesses O, proposes the partial split sigma|O (non-trivial on O).
Checked, exactly (all innovation outcomes x all legibility patterns):
  (a) every proposal on O is the restriction of a complete-data proposal (classes only shrink);
  (b) lower bound  mt_sigma >= (1 - 4 delta) m_sigma  for every sigma;
  (c) upper bound  mt_tau <= m_tau + 2 sum_{s != tau} m_s delta^{d(s,tau)}  for every tau,
      d(s, tau) = min(|A_s symdiff A_tau|, |A_s symdiff B_tau|);
  (d) when the sufficient condition of the corollary holds, the true split set is the unique maximum of mt
      (brute force over compatible sets). Exits non-zero on failure. Usage: fragmentary_check.py OUT.json"""
import sys, json, itertools, random, math

random.seed(5)


def rooted_tree(n):
    nodes = [frozenset({i}) for i in range(n)]; cl = []
    while len(nodes) > 2:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; cl.append(c); nodes.append(c)
    return [frozenset({i}) for i in range(n)] + cl + nodes


def classes(n, clades, out):
    order = sorted(range(len(clades)), key=lambda k: -len(clades[k])); val = {}
    for leaf in range(n):
        v = 'o'
        for k in order:
            if leaf in clades[k] and out[k] != 'none': v = 'a' if out[k] == 'a' else ('f', k)
        val[leaf] = v
    C = {}
    for l, v in val.items(): C.setdefault(v, set()).add(l)
    return [frozenset(c) for c in C.values()]


def compatible(A, B, X):
    Ac, Bc = X - A, X - B
    return not (A & B and A & Bc and Ac & B and Ac & Bc)


def all_splits(n):
    X = frozenset(range(n))
    return [frozenset(A) for r in range(2, n - 1) for A in itertools.combinations(range(1, n), r)]


def compatible_sets(C, X):
    out = [frozenset()]
    def rec(i, cur):
        for j in range(i, len(C)):
            if all(compatible(C[j], t, X) for t in cur):
                nxt = cur + [C[j]]; out.append(frozenset(nxt)); rec(j + 1, nxt)
    rec(0, []); return out


ok = True; rep = []; cond_cases = 0
for trial in range(40):
    n = 5; X = frozenset(range(n)); clades = rooted_tree(n); E = len(clades)
    canon = lambda A: frozenset(A) if 0 not in A else X - frozenset(A)
    true = {canon(c) for c in clades if 2 <= len(c) <= n - 2}
    rates = [random.uniform(0.01, 0.25) for _ in clades]; al = random.random()
    delta = random.uniform(0.0, 0.5); dl = [random.uniform(0, delta) for _ in range(n)]
    S = all_splits(n)
    m = {s: 0.0 for s in S}; mt = {s: 0.0 for s in S}
    pats = []
    for O in itertools.product([0, 1], repeat=n):
        pr = math.prod((1 - dl[i]) if O[i] else dl[i] for i in range(n)); pats.append((frozenset(i for i in range(n) if O[i]), pr))
    for out in itertools.product(['none', 'a', 'f'], repeat=E):
        pr = math.prod((1 - rates[k]) if o == 'none' else rates[k] * (al if o == 'a' else 1 - al) for k, o in enumerate(out))
        if pr < 1e-14: continue
        cls = classes(n, clades, out)
        full = {canon(c) for c in cls if 2 <= len(c) <= n - 2}
        for s in full: m[s] += pr
        for O, po in pats:
            props = set()
            for c in cls:
                A, B = c & O, O - c
                if len(A) >= 2 and len(B) >= 2: props.add(frozenset(A) if 0 not in A else frozenset(B))   # partial split, side w/o 0 if 0 in O
            # (a) each proposal is the restriction of a complete-data proposal
            restr = set()
            for s in full:
                A, B = s & O, (X - s) & O
                if len(A) >= 2 and len(B) >= 2: restr.add((A, B))
            pp = {(c & O, O - c) for c in cls if len(c & O) >= 2 and len(O - c) >= 2}
            ok &= all((a, b) in restr or (b, a) in restr for a, b in pp)
            # restriction-weighted support of every full split
            for s in S:
                A, B = s & O, (X - s) & O
                if len(A) >= 2 and len(B) >= 2 and ((A, B) in pp or (B, A) in pp): mt[s] += pr * po
    d = lambda s, t: min(len(s ^ t), len(s ^ (X - t)))
    for s in S:
        ok &= mt[s] >= (1 - 4 * delta) * m[s] - 1e-12
        ok &= mt[s] <= m[s] + 2 * sum(m[x] * delta ** d(x, s) for x in S if x != s) + 1e-12
    # (d) sufficient condition => unique maximum
    cond = all(m[t] + 2 * sum(m[x] * delta ** d(x, t) for x in S if x != t) < (1 - 4 * delta) * m[s]
               for t in S if t not in true for s in true if not compatible(t, s, X))
    if cond:
        cond_cases += 1
        vals = sorted(((sum(mt[x] for x in C), C) for C in compatible_sets(S, X)), key=lambda z: -z[0])
        ok &= vals[0][1] == frozenset(true) and vals[1][0] < vals[0][0]
    rep.append(dict(n=n, delta=round(delta, 3), alpha=round(al, 3), condition_held=cond))
print('ALL PASS' if ok else 'FAIL', len(rep), 'instances;', cond_cases, 'with the sufficient condition')
json.dump(dict(passed=bool(ok), instances=rep), open(sys.argv[1], 'w'), indent=1)
sys.exit(0 if ok else 1)
