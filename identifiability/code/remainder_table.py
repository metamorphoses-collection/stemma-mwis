r"""Isolated-pair C/R with a numerical remainder (answer to review finding F2, 22 Nov 2026).

For each tradition with a usable rooted tree:
  Heinrichi  - the true stemma (sampled ancestors attached as leaves), majority reading as reference;
  Donatus    - Licciardello's stemma (groups of build_matrix.py; root a polytomy), printed lemma ('0') as reference;
  Lai        - estimated tree ((A,B),(C,G),((D,E),F)) rooted on the AB edge, Bedier's printed text as reference.
(Lawrence, Getica, Iwein have no editor's or true tree in code and are not treated.)

1. Counts, Prop. prop:alpha. Pair set P = pairs isolated in the full rooted tree. A unit counts for (i, j) only if
   i, j are legible and the pair is isolated in the tree induced on the witnesses legible there (the correct guard).
     R_ij : {i,j} is a reading class, or the other legible witnesses (>= 2) form one class without i, j;
     C_ij : {i,j} is a reading class departing from the reference.
   C/R estimates alpha^2 (Prop. prop:alpha); C/(R-C) estimates alpha^2/(1-alpha^2).
2. Remainder.
   q_e  : per-edge rate, the share of units (with the clade legible and >= 2 legible witnesses outside it) at which
          the legible members of clade C_e are exactly a departing reading class (for pendant edges: the singular
          departure rate; the w(D)/L of Prop. prop:polytomynull). These are rates per RECORDED unit.
   s_cond = sum q_e. In the small-rate limit a recorded unit carries one innovation and s_cond -> 1.
   Unconditional rates: pi_e = q_e V with V = 1 - prod(1 - q_e V) (V = Pr(a place varies)); s = sum pi_e.
   P3   : Pr(|M| >= 3) under independent edges with rates pi_e (exact Poisson-binomial), and conditional on
          variation, P3/V; bound s^3/6.
   Prop. prop:alpha relative bound B = |P| Pr(|M|>=3) / sum_P R_ij (invariant under conditioning on variation),
   evaluated (a) with the plug-in P3_cond and (b) with the empirical lower bound on P3_cond, the share of recorded
   units with >= 3 distinct departing readings (which needs >= 3 innovations).
   Actual remainder at the plug-in rates: Monte Carlo of the attraction model on the tree (complete legibility),
   with alpha^2 = the observed ratio (capped at 1), giving E sum C / E sum R - alpha^2 over P.
Usage: remainder_table.py OUT.json"""
import sys, os, json, itertools, math
from collections import Counter
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
rng = np.random.default_rng(7)
HOME = os.path.expanduser('~')


# ---------------- isolation ----------------
def bad_pairs(clades, X):
    """Pairs {i,j} of X that are NOT isolated in the rooted tree with these clades (each edge once, X not a clade)."""
    cl = [c for c in set(clades) if c and c != X]
    bad, restbad = set(), set()                      # restbad: sets t with X - t a forbidden pair
    m = len(X)
    for e in cl:
        if len(e) == 2: bad.add(e)
        if len(e) == m - 2: restbad.add(e)
    for e in cl:
        for f in cl:
            if f < e:
                d = e - f
                if len(d) == 2: bad.add(d)
                if len(d) == m - 2: restbad.add(d)
    for e, f in itertools.combinations(cl, 2):
        if not (e & f) and len(e) + len(f) == m - 2: restbad.add(e | f)
    for t in restbad: bad.add(X - t)
    return bad


def induced(clades, O):
    O = frozenset(O); out = set()
    for c in clades:
        d = c & O
        if d and d != O: out.add(d)
    return out


# ---------------- counts ----------------
def pair_counts(rows, refs, clades, X, P):
    n = len(X); R = Counter(); C = Counter(); elig = Counter(); R3 = Counter(); cache = {}
    for r, ref in zip(rows, refs):
        O = frozenset(k for k in range(n) if r[k] != '?')
        if len(O) < 4: continue
        if O not in cache:
            cache[O] = bad_pairs(induced(clades, O), O)
        bad = cache[O]
        cnt = Counter(r[k] for k in O); dep3 = len(set(cnt) - {ref}) >= 3
        for i, j in P:
            if i not in O or j not in O or frozenset({i, j}) in bad: continue
            elig[(i, j)] += 1
            pairclass = r[i] == r[j] and cnt[r[i]] == 2
            others = [r[k] for k in O if k != i and k != j]
            restclass = len(set(others)) == 1 and r[i] != others[0] and r[j] != others[0]
            if pairclass or restclass:
                R[(i, j)] += 1; R3[(i, j)] += dep3
            if pairclass and r[i] != ref: C[(i, j)] += 1
    return R, C, elig, R3


def edge_rates(rows, refs, clades, n):
    q = []
    for c in clades:
        num = den = 0
        for r, ref in zip(rows, refs):
            if any(r[k] == '?' for k in c): continue
            out = [k for k in range(n) if k not in c and r[k] != '?']
            if len(out) < 2: continue
            den += 1
            v = {r[k] for k in c}
            if len(v) == 1:
                v = v.pop()
                if v != ref and all(r[k] != v for k in out): num += 1
        q.append(num / den if den else 0.0)
    return np.array(q)


def solve_V(q):
    f = lambda V: 1 - np.prod(1 - np.minimum(q * V, 1)) - V
    if q.sum() <= 1: return None
    lo, hi = 1e-9, 1.0
    if f(hi) > 0: return 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2


def p_atleast3(pi):
    d = np.zeros(4); d[0] = 1.0                     # distribution of min(|M|,3)
    for p in pi:
        nd = d * (1 - p); nd[1:] += d[:-1] * p; nd[3] += d[3] * p; d = nd
    return float(d[3])


def simulate(clades, n, pi, a2, P, N=400_000, batch=200_000):
    al = math.sqrt(min(a2, 1.0)); order = sorted(range(len(clades)), key=lambda k: -len(clades[k]))
    SR = SC = SR3 = 0.0; S3 = Sv = 0
    for b in range(0, N, batch):
        m = min(batch, N - b)
        inn = rng.random((m, len(clades))) < pi[None, :]
        nM = inn.sum(1); Sv += int((nM >= 1).sum())
        keep = nM >= 2
        inn = inn[keep]; k = inn.shape[0]; m3 = nM[keep] >= 3
        attr = rng.random((k, len(clades))) < al
        val = np.zeros((k, n), dtype=np.int32)
        for e in order:
            code = np.where(attr[:, e], 1, 2 + e)
            idx = list(clades[e])
            val[:, idx] = np.where(inn[:, e][:, None], code[:, None], val[:, idx])
        dep = [len(set(v) - {0}) for v in val[:min(k, 50000)]]
        S3 += sum(d >= 3 for d in dep) * (k / max(1, min(k, 50000)))
        for i, j in P:
            vi, vj = val[:, i], val[:, j]
            same = (val == vi[:, None]).sum(1)
            pc = (vi == vj) & (same == 2)
            oth = np.delete(val, [i, j], axis=1)
            rc = (oth == oth[:, :1]).all(1) & (vi != oth[:, 0]) & (vj != oth[:, 0])
            SR += float((pc | rc).sum()); SC += float((pc & (vi != 0)).sum()); SR3 += float(((pc | rc) & m3).sum())
    return SC / SR if SR else None, S3 / Sv if Sv else None, SR3 / SR if SR else None


# ---------------- data ----------------
def heinrichi():
    H = os.path.join(HOME, 'casc_heinrichi'); Lc = open(os.path.join(H, 'correct.txt')).read().split('\n')
    nn = int(Lc[0]); labels = Lc[1:1 + nn]; M = [l.split() for l in Lc[1 + nn:1 + 2 * nn]]
    children = {a: [labels[j] for j in range(nn) if M[i][j] == '1'] for i, a in enumerate(labels)}
    parent = {c: a for a, cs in children.items() for c in cs}; root = [a for a in labels if a not in parent][0]
    wl = [l.split() for l in open(os.path.join(H, 'heinrichi.num')).read().split('\n') if l.strip()]
    wit, rows = wl[0], [r for r in wl[1:] if len(r) == len(wl[0])]
    keep = [k for k, w in enumerate(wit) if w in labels]
    wit = [wit[k] for k in keep]; rows = [[r[k] for k in keep] for r in rows]; idx = {w: k for k, w in enumerate(wit)}
    def below(a):
        s = {idx[a]} if a in idx else set()
        for c in children[a]: s |= below(c)
        return s
    cl = set()
    for a in labels:
        if a == root: continue
        s = frozenset(below(a))
        if s: cl.add(s)
        if a in idx and children[a]: cl.add(frozenset({idx[a]}))
    for k in range(len(wit)): cl.add(frozenset({k}))
    X = frozenset(range(len(wit))); cl.discard(X)
    refs = [(Counter(v for v in r if v != '?').most_common(1) or [(None, 0)])[0][0] for r in rows]
    return wit, rows, refs, sorted(cl, key=len)


def donatus():
    sys.path.insert(0, os.path.join(HOME, 'code/stemma-mwis/code/donatus'))
    from build_matrix import GROUPS
    wl = [l.split() for l in open(os.path.join(HOME, 'stemma_xfer/donatus_final.num')).read().split('\n') if l.strip()]
    wit, rows = wl[0], [r for r in wl[1:] if len(r) == len(wl[0])]; ix = {w: k for k, w in enumerate(wit)}
    cl = {frozenset({k}) for k in range(len(wit))} | {frozenset(ix[w] for w in ms) for ms in GROUPS.values()}
    return wit, rows, ['0'] * len(rows), sorted(cl, key=len)


def lai():
    W = 'ABCDEFG'; ix = {c: k for k, c in enumerate(W)}; X = frozenset(range(7))
    U = json.load(open(os.path.join(HOME, 'MC/projects/priscus_excerpta/data/historic/lai_ombre/units_v2.json')))
    rows, refs = [], []
    for u in U:
        r = ['0'] * 7
        for k, g in enumerate(u['groups']):
            for c in g:
                if c in W: r[ix[c]] = str(k)
        rows.append(r); refs.append(str(u['lemma_group']))
    AB = frozenset({0, 1})
    sides = [frozenset({k}) for k in range(7)] + [frozenset(ix[c] for c in s) for s in ('AB', 'CG', 'DE', 'DEF')]
    cl = set()
    for s in sides:
        for t in (s, X - s):
            if t and t != X and (t <= AB or t <= X - AB): cl.add(t)
    return list(W), rows, refs, sorted(cl, key=len)


out = {}
for name, loader in (('Heinrichi', heinrichi), ('Donatus', donatus), ("Lai (AB-edge root)", lai)):
    wit, rows, refs, clades = loader(); n = len(wit); X = frozenset(range(n))
    bad = bad_pairs(clades, X)
    P = [(i, j) for i, j in itertools.combinations(range(n), 2) if frozenset({i, j}) not in bad]
    R, C, elig, R3 = pair_counts(rows, refs, clades, X, P)
    SR, SCt = sum(R.values()), sum(C.values())
    Pe = [p for p in P if elig[p]]
    rate_R = sum(R[p] / elig[p] for p in Pe)
    ratio = SCt / SR if SR else None
    q = edge_rates(rows, refs, clades, n); s_cond = float(q.sum())
    V = solve_V(q)
    pi = q * V if V else q.copy()        # no root (s_cond <= 1): V undetermined; pi = q is the largest admissible (V = 1)
    s = float(pi.sum()); P3 = p_atleast3(pi); Vm = 1 - float(np.prod(1 - pi)); P3c = P3 / Vm
    dep3 = []; rd3 = []
    for r, ref in zip(rows, refs):
        vals = [v for v in r if v != '?']
        dep3.append(len(set(vals) - {ref}) >= 3); rd3.append(len(set(vals)) >= 3)
    emp3 = sum(dep3) / len(rows); emp_r3 = sum(rd3) / len(rows)
    B_plug = len(Pe) * P3c / rate_R if rate_R else None
    B_emp = len(Pe) * emp3 / rate_R if rate_R else None
    a2 = ratio if ratio is not None else 0.5
    sim_ratio, sim_dep3, sim_f3 = simulate(clades, n, pi, a2, P if n < 20 else [p for p in P if elig[p]])
    d = dict(witnesses=n, units=len(rows), edges=len(clades), isolated_pairs_full_tree=len(P),
             pairs_with_eligible_units=len(Pe), eligible_pair_units=sum(elig.values()),
             R=SR, C=SCt, C_over_R=round(ratio, 4) if ratio is not None else None,
             C_over_R_minus_C=round(SCt / (SR - SCt), 4) if SR > SCt else None,
             pairs_with_events=sum(1 for p in P if R[p]), sum_R_rates=round(rate_R, 5),
             s_cond=round(s_cond, 4), V_hat=round(V, 4) if V else None, s_uncond=round(s, 4),
             s3_over_6=round(s ** 3 / 6, 5), P3_uncond=round(P3, 5), P3_cond_plugin=round(P3c, 4),
             share_units_3plus_readings=round(emp_r3, 3), share_units_3plus_departing=round(emp3, 3),
             rel_bound_plugin=round(B_plug, 2) if B_plug is not None else None,
             rel_bound_empirical=round(B_emp, 2) if B_emp is not None else None,
             sim_alpha2=round(min(a2, 1.0), 4), sim_ratio=round(sim_ratio, 4) if sim_ratio else None,
             sim_bias=round(sim_ratio - min(a2, 1.0), 4) if sim_ratio else None,
             R_events_at_3plus_departing=sum(R3.values()), f3_empirical_lower=round(sum(R3.values()) / SR, 4) if SR else None,
             sim_f3_R_events_with_3plus_innovations=round(sim_f3, 4) if sim_f3 is not None else None,
             sim_share_3plus_departing_given_variation=round(sim_dep3, 3) if sim_dep3 is not None else None,
             top_pairs=sorted(([wit[i] + '-' + wit[j], R[(i, j)], C[(i, j)], elig[(i, j)]] for i, j in P if R[(i, j)]),
                              key=lambda x: -x[1])[:10])
    out[name] = d
    print(name, json.dumps({k: v for k, v in d.items() if k != 'top_pairs'}))
    print('   top pairs (pair, R, C, eligible units):', d['top_pairs'])
json.dump(out, open(sys.argv[1], 'w'), indent=1)
