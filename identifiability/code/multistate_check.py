r"""Machine check of the proposition on several attractive readings (prop:multistate).
On random rooted binary trees (n = 5, 6), each edge innovates with probability pi_e; an innovation takes attractive
reading j with probability a_j (k = 2 or 3 readings) or a reading of its own. Exact expected supports of every split
(unrooted coding) and every cluster (rooted coding) are computed by enumerating all outcomes, and compared with the
attraction model at alpha = sqrt(A), A = sum a_j^2:
  * the parts with at most two innovations agree exactly (to 1e-12);
  * the full supports differ by at most P(three or more innovations);
  * the rooted support of each side of a non-trivial split is at most the unrooted support of the split.
Exits non-zero on failure. Usage: multistate_check.py OUT.json"""
import sys, json, itertools, random, math
random.seed(23)


def rooted_tree(n):
    nodes = [frozenset({i}) for i in range(n)]; cl = []
    while len(nodes) > 2:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; cl.append(c); nodes.append(c)
    return [frozenset({i}) for i in range(n)] + cl


def proposals(n, clades, outcome):
    """outcome[k]: None (no innovation), ('a', j) attractive j, ('f',) own reading."""
    X = frozenset(range(n)); order = sorted(range(len(clades)), key=lambda k: -len(clades[k])); cls = {}
    for leaf in range(n):
        v = ('o',)
        for k in order:
            if leaf in clades[k] and outcome[k] is not None:
                v = outcome[k] if outcome[k][0] == 'a' else ('f', k)
        cls.setdefault(v, set()).add(leaf)
    U, R = set(), set()
    for v, C in cls.items():
        C = frozenset(C)
        if 2 <= len(C) <= n - 2: U.add(C if 0 not in C else X - C)
        if v != ('o',) and 2 <= len(C) <= n - 1: R.add(C)
    return U, R


def supports(n, clades, rates, probs):
    """probs: list of (value, probability) for an innovation's reading. Returns dicts of exact supports, split by
    number of innovations (<=2 vs >=3)."""
    E = len(clades); lo = {}, {}; hi = {}, {}; p3 = 0.0
    choices = [None] + [v for v, _ in probs]
    pv = dict(probs)
    for out in itertools.product(choices, repeat=E):
        w = 1.0; m = 0
        for k, o in enumerate(out):
            if o is None: w *= 1 - rates[k]
            else: w *= rates[k] * pv[o]; m += 1
        if w == 0: continue
        U, R = proposals(n, clades, out)
        tgt = lo if m <= 2 else hi
        if m >= 3: p3 += w
        for C in U: tgt[0][C] = tgt[0].get(C, 0) + w
        for C in R: tgt[1][C] = tgt[1].get(C, 0) + w
    return lo, hi, p3


ok = True; rep = []
for trial in range(24):
    n = random.choice([5, 6]); clades = rooted_tree(n); X = frozenset(range(n))
    rates = [random.uniform(0.02, 0.25) for _ in clades]
    k = random.choice([2, 3]); raw = [random.random() for _ in range(k + 1)]; s = sum(raw)
    a = [x / s for x in raw[:k]]; own = raw[k] / s; A = sum(x * x for x in a); al = math.sqrt(A)
    multi = [(('a', j), a[j]) for j in range(k)] + [(('f',), own)]
    attr = [(('a', 0), al), (('f',), 1 - al)]
    loM, hiM, p3M = supports(n, clades, rates, multi)
    loA, hiA, p3A = supports(n, clades, rates, attr)
    worst_lo = worst_tr = 0.0; dom_ok = True
    for c in (0, 1):
        for C in set(loM[c]) | set(loA[c]):
            worst_lo = max(worst_lo, abs(loM[c].get(C, 0) - loA[c].get(C, 0)))
        for C in set(loM[c]) | set(loA[c]) | set(hiM[c]) | set(hiA[c]):
            mM = loM[c].get(C, 0) + hiM[c].get(C, 0); mA = loA[c].get(C, 0) + hiA[c].get(C, 0)
            worst_tr = max(worst_tr, abs(mM - mA) - p3A)
    for C in set(loM[1]) | set(hiM[1]):
        if 2 <= len(C) <= n - 2:
            canon = C if 0 not in C else X - C
            mr = loM[1].get(C, 0) + hiM[1].get(C, 0); mu = loM[0].get(canon, 0) + hiM[0].get(canon, 0)
            dom_ok &= mr <= mu + 1e-12
    t_ok = worst_lo < 1e-12 and worst_tr <= 1e-10 and abs(p3M - p3A) < 1e-10 and dom_ok
    ok &= t_ok
    rep.append(dict(n=n, k=k, A=round(A, 4), two_innovation_max_diff=worst_lo, transfer_excess=worst_tr, dominance=dom_ok, ok=t_ok))
print(f"{len(rep)} trees: two-innovation parts equal, |m - m(sqrt A)| <= P3, dominance: {'pass' if ok else 'FAIL'}")
json.dump(dict(ok=ok, trials=rep), open(sys.argv[1], 'w'), indent=1)
sys.exit(0 if ok else 1)
