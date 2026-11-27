r"""Machine check of the fragmentary-witness theorem. Complete-data units follow the attraction model on a random
rooted binary tree (n = 5); each witness is then illegible independently with probability delta_i <= delta <= 1/2,
independently of the readings. The restriction-weighted support of a full split sigma is the probability that some
reading class, restricted to the legible witnesses O, proposes the partial split sigma|O (non-trivial on O).
Checked, exactly (all innovation outcomes x all legibility patterns):
  (a) every proposal on O is the restriction of a complete-data proposal (classes only shrink);
  (b) lower bound  mt_sigma >= (1 - 4 delta) m_sigma  for every sigma;
  (c) upper bound  mt_tau <= m_tau + 2 sum_{s != tau} m_s delta^{d(s,tau)}  for every tau,
      d(s, tau) = min(|A_s symdiff A_tau|, |A_s symdiff B_tau|);
  (d) when the refined sufficient condition holds, the true split set is the unique maximum of mt
      (brute force over compatible sets);
  (e) the printed scalar bound of (c): with rho = min m_s - m_t over conflicting (false t, true s) and mu = sum m,
      delta < rho/(4+2mu) implies the refined condition; and every instance with rho > 0 is re-run at
      delta = 0.9 rho/(4+2mu), where (a)-(b) must hold and the true set must be the unique maximum of mt.
No outcome is skipped (exact enumeration). Exits non-zero on failure. Usage: fragmentary_check.py OUT.json"""
import sys, json, itertools, random, math

random.seed(5)


def rooted_tree(n):
    nodes = [frozenset({i}) for i in range(n)]; cl = []
    while len(nodes) > 2:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; cl.append(c); nodes.append(c)
    # each edge once: n pendant edges + the n-2 internal clades (the root's two child edges are among them)
    clades = [frozenset({i}) for i in range(n)] + cl
    assert len(clades) == 2 * n - 2 and len(set(clades)) == len(clades), clades
    return clades


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


def supports(n, clades, rates, al, dl, S):
    """Exact m (complete data) and mt (restriction-weighted), all innovation outcomes x all legibility patterns.
    Also returns whether (a) held. No outcome is skipped."""
    X = frozenset(range(n)); canon = lambda A: frozenset(A) if 0 not in A else X - frozenset(A)
    m = {s: 0.0 for s in S}; mt = {s: 0.0 for s in S}; a_ok = True
    pats = []
    for O in itertools.product([0, 1], repeat=n):
        pr = math.prod((1 - dl[i]) if O[i] else dl[i] for i in range(n)); pats.append((frozenset(i for i in range(n) if O[i]), pr))
    for out in itertools.product(['none', 'a', 'f'], repeat=len(clades)):
        pr = math.prod((1 - rates[k]) if o == 'none' else rates[k] * (al if o == 'a' else 1 - al) for k, o in enumerate(out))
        cls = classes(n, clades, out)
        full = {canon(c) for c in cls if 2 <= len(c) <= n - 2}
        for s in full: m[s] += pr
        for O, po in pats:
            # (a) each proposal is the restriction of a complete-data proposal
            restr = set()
            for s in full:
                A, B = s & O, (X - s) & O
                if len(A) >= 2 and len(B) >= 2: restr.add((A, B))
            pp = {(c & O, O - c) for c in cls if len(c & O) >= 2 and len(O - c) >= 2}
            a_ok &= all((a, b) in restr or (b, a) in restr for a, b in pp)
            # restriction-weighted support of every full split
            for s in S:
                A, B = s & O, (X - s) & O
                if len(A) >= 2 and len(B) >= 2 and ((A, B) in pp or (B, A) in pp): mt[s] += pr * po
    return m, mt, a_ok


def unique_max(mt, S, X, true):
    vals = sorted(((sum(mt[x] for x in C), C) for C in compatible_sets(S, X)), key=lambda z: -z[0])
    return vals[0][1] == frozenset(true) and vals[1][0] < vals[0][0]


ok = True; rep = []; cond_cases = 0; scalar_cases = 0; scalar_forced = 0
for trial in range(40):
    n = 5; X = frozenset(range(n)); clades = rooted_tree(n); E = len(clades)
    canon = lambda A: frozenset(A) if 0 not in A else X - frozenset(A)
    true = {canon(c) for c in clades if 2 <= len(c) <= n - 2}
    rates = [random.uniform(0.01, 0.25) for _ in clades]; al = random.random()
    delta = random.uniform(0.0, 0.5); dl = [random.uniform(0, delta) for _ in range(n)]
    S = all_splits(n)
    m, mt, a_ok = supports(n, clades, rates, al, dl, S); ok &= a_ok
    d = lambda s, t: min(len(s ^ t), len(s ^ (X - t)))
    def bounds_hold(m, mt, delta):
        r = True
        for s in S:
            r &= mt[s] >= (1 - 4 * delta) * m[s] - 1e-12
            r &= mt[s] <= m[s] + 2 * sum(m[x] * delta ** d(x, s) for x in S if x != s) + 1e-12
        return r
    ok &= bounds_hold(m, mt, delta)
    conflicts = [(t, s) for t in S if t not in true for s in true if not compatible(t, s, X)]
    # refined sufficient condition (from (b))
    cond = all(m[t] + 2 * sum(m[x] * delta ** d(x, t) for x in S if x != t) < (1 - 4 * delta) * m[s] for t, s in conflicts)
    # printed scalar bound of (c): local margin rho = min m_s - m_t over conflicting (false t, true s); mu = sum m
    rho = min(m[s] - m[t] for t, s in conflicts); mu = sum(m.values()); bound = rho / (4 + 2 * mu)
    scalar = rho > 0 and delta < bound
    if scalar:
        scalar_cases += 1; ok &= cond                 # the scalar bound implies the refined condition
    if cond:
        cond_cases += 1; ok &= unique_max(mt, S, X, true)
    # the scalar bound, forced: when rho > 0, re-run at an illegibility rate just inside rho/(4+2mu)
    forced = None
    if rho > 0:
        df = 0.9 * bound; dlf = [random.uniform(0, df) for _ in range(n)]
        _, mtf, a_okf = supports(n, clades, rates, al, dlf, S)
        forced = bool(a_okf and bounds_hold(m, mtf, df) and unique_max(mtf, S, X, true))
        ok &= forced; scalar_forced += 1
    rep.append(dict(n=n, delta=round(delta, 3), alpha=round(al, 3), rho=round(rho, 5), mu=round(mu, 4),
                    scalar_bound=round(bound, 5), scalar_held=scalar, condition_held=cond, forced_scalar_pass=forced))
print('ALL PASS' if ok else 'FAIL', len(rep), 'instances;', cond_cases, 'with the refined sufficient condition;',
      scalar_cases, 'with the printed scalar bound delta < rho/(4+2mu);', scalar_forced,
      'instances with rho > 0 re-run at delta = 0.9 rho/(4+2mu): bounds hold and the true set is the unique maximum')
json.dump(dict(passed=bool(ok), instances=rep, refined_condition_cases=cond_cases, scalar_bound_cases=scalar_cases,
               scalar_forced_cases=scalar_forced), open(sys.argv[1], 'w'), indent=1)
sys.exit(0 if ok else 1)
