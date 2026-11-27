r"""Theorem cuttingconsistency (b), checked by simulation (check C26).

The cutting-plane estimator selects a binary tree T_hat maximising W(T') = number of (unit, distinct partial
proposal) pairs displayed by T'. Part (b): with i.i.d. units, T_hat = T with probability tending to one if the true T
is the unique maximum of M(T') = E[W(T')]/L over binary trees, and with probability tending to zero if some T' has
M(T') > M(T).

Model: attraction model on a random rooted binary tree with n = 5 leaves, unrooted coding, each witness illegible
independently with probability delta (as in cutting_consistency.py, check C25). M is computed exactly by enumeration;
units are then sampled and the selection frequency of T is measured at increasing L.
Case A: instances where T is the unique maximum of M. Pass: frequency of T_hat = T at the largest L is >= 0.95, and
        larger than at the smallest L.
Case B: instances where some T' has M(T') > M(T) (long branches). Pass: frequency at the largest L is <= 0.05, and
        smaller than at the smallest L.
Ties in W are counted as failures to select T. Usage: cutting_consistency_sim.py OUT.json"""
import sys, json, random, itertools
from collections import defaultdict
random.seed(2611)
n = 5
X = frozenset(range(n))


def nontrivial_splits():
    return [frozenset(c) for r in range(2, n - 1) for c in itertools.combinations(range(1, n), r) if r <= n - 2]


def compatible(A, B):
    return not (A & B and A - B and B - A and X - (A | B))


TREES = [frozenset(c) for c in itertools.combinations(nontrivial_splits(), n - 3)
         if all(compatible(a, b) for a, b in itertools.combinations(c, 2))]


def canon(O, S):
    m = min(O); return (O, frozenset(S if m in S else O - S))


def displays(t, p):
    O, S = p
    return any((A & O) == S or (O - (A & O)) == S for A in t)


def random_clades():
    nodes = [frozenset([i]) for i in range(n)]; cl = []
    while len(nodes) > 1:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); nodes.append(a | b); cl += [a, b]
    return cl


def classes(cl, st):
    read = {}
    for x in range(n):
        best = None
        for e in range(len(cl)):
            if st[e] and x in cl[e] and (best is None or len(cl[e]) < len(cl[best])): best = e
        read[x] = 0 if best is None else ('a' if st[best] == 1 else ('o', best))
    cls = defaultdict(set)
    for x, r in read.items(): cls[r].add(x)
    return [frozenset(c) for c in cls.values()]


def proposals(cls, O):
    if len(O) < 4: return set()
    out = set()
    for C in cls:
        S = C & O
        if 2 <= len(S) <= len(O) - 2: out.add(canon(O, S))
    return out


def exact_M(cl, pi, alpha, delta):
    M = [0.0] * len(TREES)
    legs = []
    for pat in itertools.product((0, 1), repeat=n):
        O = frozenset(i for i in range(n) if pat[i]); pr = 1.0
        for i in range(n): pr *= (1 - delta) if pat[i] else delta
        legs.append((O, pr))
    for st in itertools.product((0, 1, 2), repeat=len(cl)):
        pr = 1.0
        for e, s in enumerate(st): pr *= (1 - pi[e]) if s == 0 else pi[e] * (alpha if s == 1 else 1 - alpha)
        cls = classes(cl, st)
        for O, pl in legs:
            for p in proposals(cls, O):
                for i, t in enumerate(TREES):
                    if displays(t, p): M[i] += pr * pl
    return M


def sample_unit(cl, pi, alpha, delta):
    st = []
    for e in range(len(cl)):
        u = random.random()
        st.append(0 if u >= pi[e] else (1 if random.random() < alpha else 2))
    O = frozenset(i for i in range(n) if random.random() >= delta)
    return proposals(classes(cl, st), O)


def true_tree(cl):
    return frozenset(A if 0 not in A else X - A for A in cl if 2 <= len(A) <= n - 2)


def selection_freq(cl, pi, alpha, delta, iT, L, reps):
    hit = 0
    for _ in range(reps):
        W = [0] * len(TREES)
        for _ in range(L):
            for p in sample_unit(cl, pi, alpha, delta):
                for i, t in enumerate(TREES):
                    if displays(t, p): W[i] += 1
        best = max(W)
        hit += (W[iT] == best and sum(w == best for w in W) == 1)
    return hit / reps


LS = [25, 100, 400, 1600]; REPS = 100
cases = {'A': [], 'B': []}
tries = 0
while (len(cases['A']) < 3 or len(cases['B']) < 3) and tries < 400:
    tries += 1
    cl = random_clades(); T = true_tree(cl); iT = TREES.index(T)
    long_branch = random.random() < 0.5
    pi = [random.uniform(0.02, 0.15) for _ in cl]
    if long_branch:                                          # make two non-sister pendant edges long
        leaves = [e for e, c in enumerate(cl) if len(c) == 1]
        a, b = random.sample(leaves, 2); pi[a] = pi[b] = random.uniform(0.35, 0.6)
        for e, c in enumerate(cl):
            if len(c) > 1: pi[e] = random.uniform(0.005, 0.03)
    alpha = random.uniform(0, 1); delta = random.uniform(0.02, 0.2)
    M = exact_M(cl, pi, alpha, delta)
    others = max(M[i] for i in range(len(TREES)) if i != iT)
    gap = M[iT] - others
    if gap > 0.02 and len(cases['A']) < 3:
        cases['A'].append((cl, pi, alpha, delta, iT, gap))
    elif gap < -0.02 and len(cases['B']) < 3:
        cases['B'].append((cl, pi, alpha, delta, iT, gap))

rows = []; ok = True
for kind in ('A', 'B'):
    for cl, pi, alpha, delta, iT, gap in cases[kind]:
        freqs = [selection_freq(cl, pi, alpha, delta, iT, L, REPS) for L in LS]
        if kind == 'A': passed = freqs[-1] >= 0.95 and freqs[-1] > freqs[0]
        else: passed = freqs[-1] <= 0.05 and freqs[-1] < freqs[0]
        ok &= passed
        rows.append(dict(case=kind, gap_M=round(gap, 4), alpha=round(alpha, 3), delta=round(delta, 3), L=LS,
                         freq_true_selected=freqs, passed=passed))
        print(kind, 'gap', round(gap, 4), 'freq', freqs, 'pass', passed, flush=True)
found = {k: len(v) for k, v in cases.items()}
ok &= found['A'] == 3 and found['B'] == 3
out = dict(passed=bool(ok), instances_found=found, reps=REPS, rows=rows)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
print({k: out[k] for k in ('passed', 'instances_found', 'reps')})
sys.exit(0 if ok else 1)
