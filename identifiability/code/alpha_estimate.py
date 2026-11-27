r"""Estimating the coincidence probability alpha apart from descent and from the selection of variant places.

For a pair of leaves (i, j), compare in the same collation
  R_ij = units where {i, j} is a reading class or the other legible witnesses form one class without them
         (the unrooted coding proposes {i,j}|rest: retention or coincidence), and
  C_ij = units where {i, j} is a reading class that departs from the reference (the rooted coding proposes {i,j}).
The pair is *isolated* in a rooted tree if, at second order, only the pair of pendant edges of i and j produces these
proposals (Theorem mechanisms): {i,j} is not a clade, no clade C_e has C_e minus a sub-clade equal to {i,j}, and the
rest X minus {i,j} is neither a union of two disjoint clades, nor a clade minus a sub-clade, nor a clade. For isolated
pairs the two-mark parts satisfy C = alpha^2 R exactly, so alpha^2 is estimated by sum C / sum R over isolated pairs;
units that do not vary propose nothing in either coding, so the selection of variant places cancels.
Part 1 checks the identity C2 = alpha^2 R2 for the two-mark parts on random trees (exhaustive). Part 2 estimates
alpha^2 on Heinrichi (true tree, majority reading as reference) and on the Lai (estimated tree, pairs isolated under
every rooting, Bedier's printed text as reference). Usage: alpha_estimate.py OUT.json"""
import sys, os, json, itertools, random, math
from collections import Counter

random.seed(13)


def isolated(i, j, clades, X):
    pair, rest = frozenset({i, j}), X - frozenset({i, j})
    cl = set(clades)
    if pair in cl: return False
    for e in clades:
        for f in clades:
            if f < e and (e - f) in (pair, rest): return False
    for e, f in itertools.combinations(clades, 2):
        if {e, f} == {frozenset({i}), frozenset({j})}: continue      # the intended pair of pendant edges
        if not (e & f) and (e | f) in (pair, rest): return False
    return rest not in cl


# ---------- part 1: exact two-mark identity on random rooted trees ----------
def rooted_tree(n):
    nodes = [frozenset({i}) for i in range(n)]; cl = []
    while len(nodes) > 2:
        a, b = random.sample(nodes, 2); nodes.remove(a); nodes.remove(b); c = a | b; cl.append(c); nodes.append(c)
    # each edge once: n pendant edges + the n-2 internal clades (the root's two child edges are among them)
    clades = [frozenset({i}) for i in range(n)] + cl
    assert len(clades) == 2 * n - 2 and len(set(clades)) == len(clades), clades
    return clades


def proposals(n, clades, out):
    X = frozenset(range(n)); order = sorted(range(len(clades)), key=lambda k: -len(clades[k])); val = {}
    for leaf in range(n):
        v = 'o'
        for k in order:
            if leaf in clades[k] and out[k] != 'none': v = 'a' if out[k] == 'a' else ('f', k)
        val[leaf] = v
    C = {}
    for l, v in val.items(): C.setdefault(v, set()).add(l)
    return {(frozenset(c), v) for v, c in C.items()}


ok = True; tested = 0
for _ in range(60):
    n = random.choice([5, 6, 7]); X = frozenset(range(n)); clades = rooted_tree(n); E = len(clades)
    rates = [random.uniform(0.02, 0.3) for _ in clades]; al = random.random()
    for i, j in itertools.combinations(range(n), 2):
        if not isolated(i, j, clades, X): continue
        pair, rest = frozenset({i, j}), X - frozenset({i, j})
        R2 = C2 = 0.0
        for k in range(E):                                    # one mark never proposes these (they are not clades)
            for v in 'af':
                out = ['none'] * E; out[k] = v; P = proposals(n, clades, out)
                ok &= not any(c in (pair, rest) for c, _ in P)
        for e, f in itertools.combinations(range(E), 2):
            base = rates[e] * rates[f] * math.prod(1 - rates[g] for g in range(E) if g not in (e, f))
            for ve, vf in itertools.product('af', repeat=2):
                w = (al if ve == 'a' else 1 - al) * (al if vf == 'a' else 1 - al)
                out = ['none'] * E; out[e] = ve; out[f] = vf; P = proposals(n, clades, out)
                if any(c == pair or (c == rest and len(rest) >= 2) for c, _ in P): R2 += base * w
                if any(c == pair and v != 'o' for c, v in P): C2 += base * w
        ok &= abs(C2 - al * al * R2) < 1e-12; tested += 1
print('part 1:', 'PASS' if ok else 'FAIL', tested, 'isolated pairs')


# ---------- part 2: data ----------
def counts(wit, rows, ref, pairs, complete_only=False, guard=None):
    """R and C for each pair; ref(r) returns the reference reading of a unit (or None).
    complete_only: use only units in which every witness is legible. Otherwise classes are computed among the
    legible witnesses, and an inherited family reading whose other members are illegible can look like an
    agreement of an isolated pair (review of 18 Nov 2026, finding 2): isolation is a property of the full tree."""
    out = {}
    for i, j in pairs:
        R = C = 0
        for r in rows:
            if r[i] == '?' or r[j] == '?': continue
            if complete_only and '?' in r: continue
            if guard is not None and any(r[k] == '?' for k in guard[(i, j)]): continue
            O = [k for k in range(len(wit)) if r[k] != '?']
            cls = Counter(r[k] for k in O)
            pairclass = r[i] == r[j] and cls[r[i]] == 2
            others = [r[k] for k in O if k not in (i, j)]
            restclass = len(others) >= 2 and len(set(others)) == 1 and r[i] != others[0] and r[j] != others[0]
            if pairclass or restclass: R += 1
            if pairclass and r[i] != ref(r): C += 1
        out[(wit[i], wit[j])] = (R, C)
    return out


res = dict(part1_pass=bool(ok), part1_pairs=tested)
# Heinrichi: true tree (directed adjacency, sampled ancestors attached as leaves), majority reading as reference
H = os.path.expanduser('~/casc_heinrichi'); Lc = open(os.path.join(H, 'correct.txt')).read().split('\n')
nn = int(Lc[0]); labels = Lc[1:1 + nn]; M = [l.split() for l in Lc[1 + nn:1 + 2 * nn]]
children = {a: [labels[j] for j in range(nn) if M[i][j] == '1'] for i, a in enumerate(labels)}
parent = {c: a for a, cs in children.items() for c in cs}
root = [a for a in labels if a not in parent][0]
wl = [l.split() for l in open(os.path.join(H, 'heinrichi.num')).read().split('\n') if l.strip()]
wit, rows = wl[0], [r for r in wl[1:] if len(r) == len(wl[0])]
keep = [k for k, w in enumerate(wit) if w in labels]
wit = [wit[k] for k in keep]; rows = [[r[k] for k in keep] for r in rows]
idx = {w: k for k, w in enumerate(wit)}
def below(a):
    s = {idx[a]} if a in idx else set()
    for c in children[a]: s |= below(c)
    return s
Xh = frozenset(range(len(wit)))
cladesH = []
for a in labels:
    if a == root: continue
    s = frozenset(below(a))
    if s: cladesH.append(s)
    if a in idx and children[a]: cladesH.append(frozenset({idx[a]}))       # pendant edge of a sampled ancestor
cladesH = list(set(cladesH))
pairsH = [(i, j) for i, j in itertools.combinations(range(len(wit)), 2) if isolated(i, j, cladesH, Xh)]
def majority(r):
    c = Counter(v for v in r if v != '?'); return c.most_common(1)[0][0]
cH = counts(wit, rows, majority, pairsH)
RH, CH = sum(v[0] for v in cH.values()), sum(v[1] for v in cH.values())
res['heinrichi'] = dict(witnesses=len(wit), isolated_pairs=len(pairsH), R=RH, C=CH, alpha2_hat=round(CH / RH, 4) if RH else None)
cHc = counts(wit, rows, majority, pairsH, complete_only=True)
# local guard: every witness of the smallest true clade containing i and j must be legible (an inherited reading
# can pose as the pair {i,j} only through a clade containing both, whose other members are then visible)
allH = frozenset(range(len(wit)))
gH = {(i, j): min((c for c in cladesH if i in c and j in c), key=len, default=allH) for i, j in pairsH}
cHg = counts(wit, rows, majority, pairsH, guard=gH)
RHg, CHg = sum(v[0] for v in cHg.values()), sum(v[1] for v in cHg.values())
res['heinrichi_local_guard'] = dict(R=RHg, C=CHg, alpha2_hat=round(CHg / RHg, 4) if RHg else None,
                                    pairs_with_counts=sum(1 for v in cHg.values() if v[0]))
print('Heinrichi, local legibility guard:', res['heinrichi_local_guard'])
RHc, CHc = sum(v[0] for v in cHc.values()), sum(v[1] for v in cHc.values())
res['heinrichi_complete_units_only'] = dict(units=sum('?' not in r for r in rows), units_all=len(rows), R=RHc, C=CHc,
                                           alpha2_hat=round(CHc / RHc, 4) if RHc else None)
print('Heinrichi, complete units only:', res['heinrichi_complete_units_only'])
print('Heinrichi:', res['heinrichi'])
# Lai: estimated unrooted tree ((A,B),(C,G),((D,E),F)); pairs isolated under every rooting; printed text as reference
W = 'ABCDEFG'; Xl = frozenset(range(7)); ix = {c: k for k, c in enumerate(W)}
splits = [frozenset(ix[c] for c in s) for s in ('AB', 'CG', 'DE', 'DEF')]
edges_sides = [frozenset({k}) for k in range(7)] + splits
def clades_rooted_at(side):            # root on the edge with this side: clades are sides away from the root
    cl = []
    for s in edges_sides:
        for t in (s, Xl - s):
            if len(t) == 0 or len(t) == 7: continue
            # t is a clade iff it does not contain the root edge's far side: t subset of side or of complement, nested
            if (t <= side or t <= Xl - side): cl.append(t)
    return list(set(cl))
rootings = {'Bedier (root on the DEF edge)': clades_rooted_at(frozenset(ix[c] for c in 'DEF')),
            'AB against the rest': clades_rooted_at(frozenset(ix[c] for c in 'AB')),
            'Quentin (central vertex)': list({frozenset(ix[c] for c in t) for t in ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'AB', 'CG', 'DE', 'DEF')})}
pairsL = sorted({(i, j) for i, j in itertools.combinations(range(7), 2)
                 if any(isolated(i, j, cl, Xl) for cl in rootings.values())})
U = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'lai_ombre', 'units_v2.json')))
rowsL, refs = [], []
for u in U:
    r = ['0'] * 7
    for k, g in enumerate(u['groups']):
        for c in g:
            if c in W: r[ix[c]] = str(k)
    rowsL.append(r); refs.append(str(u['lemma_group']))
refmap = {id(r): f for r, f in zip(rowsL, refs)}
cL = counts(list(W), rowsL, lambda r: refmap[id(r)], pairsL)
RL, CL = sum(v[0] for v in cL.values()), sum(v[1] for v in cL.values())
per_root = {}
for name, cl in rootings.items():
    pr = [(i, j) for i, j in pairsL if isolated(i, j, cl, Xl)]
    Rr = sum(cL[(W[i], W[j])][0] for i, j in pr); Cr = sum(cL[(W[i], W[j])][1] for i, j in pr)
    per_root[name] = dict(pairs=[W[i] + W[j] for i, j in pr], R=Rr, C=Cr, alpha2_hat=round(Cr / Rr, 4) if Rr else None)
    print('Lai,', name, per_root[name])
res['lai_by_rooting'] = per_root
res['lai'] = dict(isolated_pairs=[''.join(p) for p in cL], per_pair={''.join(k): v for k, v in cL.items()}, R=RL, C=CL,
                  alpha2_hat=round(CL / RL, 4) if RL else None)
print('Lai:', res['lai'])
json.dump(res, open(sys.argv[1], 'w'), indent=1)
sys.exit(0 if ok else 1)
