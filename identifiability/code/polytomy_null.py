r"""Proposition polytomynull: the coincidence null for a pairing at a polytomy, checked exactly.

Attraction model on a rooted tree whose leaves are the witnesses: every non-root vertex has an edge with innovation
probability pi_e; an innovation takes the attractive reading with probability alpha, otherwise a reading of its own;
a witness carries the reading of the lowest innovating edge on its path (or the archetype's). Rooted coding: a unit
proposes the set of witnesses sharing each non-archetype reading. Let v have children with clades C_1..C_k, k >= 3,
and tau = C_i u C_j. Claim:
  m(tau) = alpha^2 m(C_i) m(C_j) + [k = 3, v not the root] (1 - alpha^2) m(C_v) m(C_l) + O(eps^3),
with every pi of order eps, where C_l is the third clade and m(C) is the probability that a unit proposes C.
Exact check: all 3^|E| innovation states are enumerated on random trees with a polytomy, at eps = 0.004 and 0.002 (total rate up to about 0.06); the
remainder must scale as eps^3 (ratio near 8); its size relative to the leading term is reported (it grows like eps / alpha^2). Usage: polytomy_null.py OUT.json"""
import sys, json, random, itertools
from collections import defaultdict
random.seed(11)


def random_tree_with_polytomy(k):
    """returns parent list (vertex 0 = root), leaves, the polytomy vertex and its children."""
    parent = [-1]
    def grow(p, depth):
        v = len(parent); parent.append(p)
        if depth > 0 and random.random() < 0.7:
            for _ in range(2): grow(v, depth - 1)
        return v
    # root -> optional stem -> v (polytomy) with k children subtrees; plus an outgroup subtree at the root
    root_is_v = random.random() < 0.4
    if root_is_v: v = 0
    else:
        v = len(parent); parent.append(0); grow(0, 1)          # v and an outgroup
    kids = [grow(v, random.choice([0, 0, 1])) for _ in range(k)]
    return parent, v, kids


def clades(parent):
    n = len(parent); ch = defaultdict(list)
    for x, p in enumerate(parent):
        if p >= 0: ch[p].append(x)
    leaves = [x for x in range(n) if not ch[x]]
    below = {}
    for x in reversed(range(n)):
        below[x] = frozenset([x]) if not ch[x] else frozenset().union(*(below[c] for c in ch[x]))
    return leaves, below


def exact_m(parent, pi, alpha):
    leaves, below = clades(parent)
    E = [x for x in range(len(parent)) if parent[x] >= 0]
    m = defaultdict(float)
    for st in itertools.product((0, 1, 2), repeat=len(E)):     # 0 none, 1 attractive, 2 own
        pr = 1.0
        for e, s in zip(E, st): pr *= (1 - pi[e]) if s == 0 else pi[e] * (alpha if s == 1 else 1 - alpha)
        if pr == 0: continue
        mark = dict(zip(E, st)); read = {}
        for lf in leaves:
            x, r = lf, 0
            while x > 0 or (x == 0 and False):
                if mark.get(x, 0): r = 'a' if mark[x] == 1 else x; break
                x = parent[x]
                if x == 0: break
            read[lf] = r
        groups = defaultdict(set)
        for lf, r in read.items():
            if r != 0: groups[r].add(lf)
        for g in groups.values(): m[frozenset(g)] += pr
    return m, below


res = []; worst = 0.0; ok = True
for trial in range(60):
    k = random.choice([3, 3, 4])
    parent, v, kids = random_tree_with_polytomy(k)
    if len(parent) > 12: continue
    alpha = random.uniform(0.2, 0.95); base = {x: random.uniform(0.5, 2.0) for x in range(len(parent))}
    i, j = random.sample(range(k), 2); l = [t for t in range(k) if t not in (i, j)]
    rems = []
    for eps in (0.004, 0.002):
        pi = {x: base[x] * eps for x in range(len(parent))}
        m, below = exact_m(parent, pi, alpha)
        tau = below[kids[i]] | below[kids[j]]
        lead = alpha ** 2 * m[below[kids[i]]] * m[below[kids[j]]]
        if k == 3 and v != 0: lead += (1 - alpha ** 2) * m[below[v]] * m[below[kids[l[0]]]]
        rems.append((m[tau] - lead, lead, m[tau]))
    (r1, l1, t1), (r2, l2, t2) = rems
    ratio = abs(r1 / r2) if r2 else float('inf')
    rel = abs(r1) / l1
    passed = abs(r1) < 1e-14 or 7 <= ratio <= 9          # remainder O(eps^3): halving eps divides it by ~8
    ok &= passed; worst = max(worst, rel)
    res.append(dict(edges=len(parent) - 1, k=k, v_is_root=v == 0, alpha=round(alpha, 3), m_tau=t1, leading=l1,
                    rel_remainder_eps0004=rel, remainder_ratio_halving=ratio, passed=passed))
out = dict(instances=len(res), all_pass=bool(ok), worst_relative_remainder_at_eps_0_004=worst,
           k_counts={kk: sum(r['k'] == kk for r in res) for kk in (3, 4)}, root_polytomies=sum(r['v_is_root'] for r in res), rows=res)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
print({x: out[x] for x in ('instances', 'all_pass', 'worst_relative_remainder_at_eps_0_004', 'k_counts', 'root_polytomies')})
print('ratios', sorted(round(r['remainder_ratio_halving'], 2) for r in res)[:5], '...', sorted(round(r['remainder_ratio_halving'], 2) for r in res)[-5:])
sys.exit(0 if ok else 1)
