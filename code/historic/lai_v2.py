"""Lai de l'Ombre v2 (22 Nov 2026): corrected units and every number of the paper's Lai section, recomputed.

Corrections to data/lai_ombre/units.json (which is left untouched):
 (1) Doubled page. n056 is the tissue guard in front of p. 1: the vision model read p. 1 THROUGH the tissue, faintly,
     and invented readings (e.g. v.1 'EG', v.3 'B ... A ...', 'prologue (v. 1-42)'). n058 is p. 1 itself; its read was
     checked against the page image (archive.org lelaidelombre00jeanuoft n58) and is correct ('prologue (v. 1-52)
     manque en F', consistent with F's first citation at v.53). All n056 units are dropped; n058 is kept.
 (2) Composite notes. A note that puts a witness in two reading groups covers several words. It is split into
     elementary units, one per lemma. Rule: every printed reading group (sigla + reading) goes into exactly ONE
     elementary unit, so no departure is created or lost; groups that are alternative readings of the same word
     share a unit. Each elementary unit satisfies exclusivity (asserted). Splits checked on the images where the
     overlap proposes crossing groups (v.176, p.10; v.256, p.13) or a doubled siglum (v.307, p.16): all as printed.
     v.176: base 'Quar m'eust cele fet'; BF 'C'or' (word 1), CDF 'ele' / G 'ore' (word 3): two loci, both real.
     v.256: BCG 'Il n'entendoit' (verb) and CF 'a ce dire' (end of verse): two loci, both real.

Outputs: data/lai_ombre/units_v2.json, data/lai_ombre/splits_v2.tsv, data/lai_ombre/results_v2.json.
Usage: python3 code/historic/lai_v2.py   (a few seconds)"""
import json, os, re, sys, itertools, collections, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(ROOT, 'data', 'lai_ombre')
W = 'ABCDEFG'; ALL = frozenset(W)
U = json.load(open(os.path.join(D, 'units.json')))

DROP_PAGES = {'n056'}
# note index in units.json -> list of elementary units; each unit = (list of indices into groups[1:], lemma label)
SPLITS = {
    12:  [([0], 'que ja le'), ([1], 'Mais p.'), ([2], 'le sen')],
    20:  [([0], "s'aucuns fous"), ([1], 'le langue')],
    26:  [([0], 'des bons'), ([1], 'ge di pieça')],
    35:  [([0], 'Et'), ([1], 'a foi le')],
    38:  [([0], "d'iluec avant"), ([1, 2, 3], 'mesure')],
    50:  [([0], 'On d.'), ([1], 'bien riue'), ([2], 'et (omission)')],
    80:  [([0], 'Le trovast on'), ([1], 'on plus que ne')],
    103: [([0, 2], 'Et frans et larges'), ([1], 'order of vv. 109-110')],
    108: [([0, 1], 'En ce'), ([2], 'corut')],
    110: [([0, 1], 'gr.'), ([2, 3], 'ot eu')],
    112: [([0], 'Ne li f.'), ([1, 2, 3], 'entrues')],
    154: [([0, 1], 'e. et de'), ([2], 'G omits the verse (as printed)')],
    155: [([0, 1, 2], 'mon cuer alaschier'), ([3], "Ne m'i puet")],
    159: [([0], "Quar -> C'or"), ([1, 2], 'cele')],
    160: [([0, 1], 'dolz b.'), ([2], 'entor le c.')],
    219: [([0, 1], "Il n'entendoit"), ([2, 3], 'a ce dire')],
    228: [([0, 1, 2, 3], 'Il dist'), ([4], 'je me t.')],
    258: [([0, 1, 3], 'de v.'), ([2], 'E de v., second entry (as printed)')],
    340: [([0], 'a vos'), ([1, 2], 'dous iex')],
    341: [([0], 'Moult'), ([1], 'orainz')],
    359: [([0], 'Et vous'), ([1], 'autre sens')],
    360: [([0, 1], 'folement'), ([2], 'si me p.')],
    387: [([0], 'Je me'), ([1], 'mis et met')],
    395: [([0, 1, 2, 3], 'un petit risset'), ([4], 'C: lines added after v.465')],
    452: [([0], 'Se je vous'), ([1], 'a la v.')],
    494: [([0, 1, 2, 3], 'se g.'), ([4], 'E second entry (as printed)')],
    551: [([0], 'dist il'), ([1], 'bien le cuit f.')],
    552: [([1, 2], 'Vostre volenté'), ([0, 3], "jusqu'a son")],
    613: [([0], 'fait ele'), ([1], 'a. ça')],
    660: [([0], 'trés'), ([1, 2], 'c. parmi'), ([3], 'la cuisse')],
    661: [([1, 2], 'Ne soufferoie je'), ([0, 3], 'tel anguisse')],
    673: [([0, 2], 'Elle dist'), ([1], 'donques')],
    707: [([0, 1], 'anieus'), ([2], 'Con d.')],
    752: [([0], 'En (omission)'), ([1, 2, 3], 'ert et b.')],
}


def segments(raw):
    """(sigla, reading) segments of a raw note, in print order, parsed as lai_extract.py does."""
    t = re.sub(r'\([^()]*\)', ' ', raw)
    m = re.match(r'^\s*(\d{1,3})(?:\s*-\s*(\d{1,3}))?\s+(.*)$', t.strip())
    body = m.group(3)
    hits = [g for g in re.finditer(r'(?:^|,\s+|;\s+)(' + r'[A-G]{1,7}' + r')(?=\s|$)', body) if len(set(g.group(1))) == len(g.group(1))]
    out = []
    for k, g in enumerate(hits):
        end = hits[k + 1].start() if k + 1 < len(hits) else len(body)
        out.append((g.group(1), body[g.end(1):end].strip(' ,;')))
    return out


def overlaps(groups):
    c = collections.Counter(w for g in groups if g for w in set(g))
    return any(v > 1 for v in c.values())


def build():
    over = [i for i, u in enumerate(U) if overlaps(u['groups'][1:])]
    assert len(over) == 35, len(over)
    assert set(over) - {i for i, u in enumerate(U) if u['page'] in DROP_PAGES} == set(SPLITS), sorted(set(over) ^ set(SPLITS))
    out, rows = [], []
    for i, u in enumerate(U):
        if u['page'] in DROP_PAGES: continue
        if i not in SPLITS:
            out.append(dict(u, src=i, sub=0, lemma_label=None)); continue
        gs = u['groups'][1:]; seg = segments(u['readings'][0])
        assert [''.join(sorted(s)) for s, _ in seg] == gs, (i, seg, gs)
        used = sorted(k for ks, _ in SPLITS[i] for k in ks)
        assert used == list(range(len(gs))), (i, used)          # every printed group in exactly one unit
        for j, (ks, lab) in enumerate(SPLITS[i]):
            g2 = [gs[k] for k in ks]
            assert not overlaps(g2), (i, j, g2)                  # exclusivity
            rd = [f'{seg[k][0]} {seg[k][1]}' for k in ks]
            out.append(dict(page=u['page'], line=u['line'], entry_line=u['entry_line'], leaf=None, lemma=None,
                            lemma_group=0, groups=[''] + g2, readings=rd, src=i, sub=j + 1, lemma_label=lab))
            rows.append((u['line'], u['page'], i, j + 1, lab, ' '.join(g2), ' | '.join(rd), u['readings'][0]))
    assert not any(overlaps(u['groups'][1:]) for u in out)
    json.dump(out, open(os.path.join(D, 'units_v2.json'), 'w'), ensure_ascii=False, indent=0)
    with open(os.path.join(D, 'splits_v2.tsv'), 'w') as f:
        f.write('verse\tpage\tnote_index_v1\tunit\tlemma\tgroups\treadings\tprinted_note_v1\n')
        for r in rows: f.write('\t'.join(str(x).replace('\t', ' ') for x in r) + '\n')
    return out, len(rows)


# ---------------------------------------------------------------- scoring
def k(f): return ''.join(sorted(f))


def unit_groups(u):
    return {frozenset(g) for g in u['groups'][1:] if g and len(g) >= 2 and frozenset(g) != ALL}


def unrooted_props(u):
    s = set()
    for g in unit_groups(u):
        if 2 <= len(g) <= 5: s.add(g if 'A' not in g else ALL - g)
    return s


def trees(leaves):
    leaves = tuple(leaves)
    if len(leaves) == 1: yield frozenset(); return
    first, rest = leaves[0], leaves[1:]
    for r in range(len(rest)):
        for comb in itertools.combinations(rest, r):
            L = (first,) + comb; R = tuple(x for x in rest if x not in comb)
            for tl in trees(L):
                for tr in trees(R):
                    cl = set(tl) | set(tr)
                    if len(L) >= 2: cl.add(frozenset(L))
                    if len(R) >= 2: cl.add(frozenset(R))
                    yield frozenset(cl)


SPL = sorted({frozenset(c) for r in (2, 3, 4, 5) for c in itertools.combinations('BCDEFG', r)}, key=lambda s: (len(s), k(s)))
SPL = [s for s in SPL if len(s) <= 5]                                      # unrooted splits keyed by the side without A
CLU = sorted({frozenset(c) for r in range(2, 7) for c in itertools.combinations(W, r)}, key=lambda s: (len(s), k(s)))
UT = [frozenset(c for c in t if 2 <= len(c) <= 5) for t in trees('BCDEFG')]   # 945 unrooted binary trees
RT = list(trees(W))                                                        # 10,395 rooted binary trees
si = {s: j for j, s in enumerate(SPL)}; ci = {c: j for j, c in enumerate(CLU)}
MU = np.zeros((len(UT), len(SPL)), np.int8)
for a, t in enumerate(UT):
    for s in t: MU[a, si[s]] = 1
MR = np.zeros((len(RT), len(CLU)), np.int8)
for a, t in enumerate(RT):
    for c in t: MR[a, ci[c]] = 1


def crosses(a, b):
    a2, b2 = ALL - a, ALL - b
    return all([a & b, a & b2, a2 & b, a2 & b2])


def unit_matrices(units):
    XU = np.zeros((len(units), len(SPL)), np.int32); XR = np.zeros((len(units), len(CLU)), np.int32)
    for r, u in enumerate(units):
        for s in unrooted_props(u): XU[r, si[s]] = 1
        for c in unit_groups(u): XR[r, ci[c]] = 1
    return XU, XR


def trinomial_p(x):
    """Prop. 7.3 extended to a trio: under equal rates, (x1,x2,x3) ~ Mult(N, 1/3); statistic = one-sided binomial
    p-value of the largest count against the second, given their sum; p = P(statistic <= observed)."""
    def tail(a, n): return sum(math.comb(n, i) for i in range(a, n + 1)) / 2 ** n
    N = sum(x); a, b, _ = sorted(x, reverse=True); s_obs = tail(a, a + b); p = 0.0
    for x1 in range(N + 1):
        for x2 in range(N - x1 + 1):
            x3 = N - x1 - x2; q = sorted((x1, x2, x3), reverse=True)
            if tail(q[0], q[0] + q[1]) <= s_obs + 1e-15:
                p += math.factorial(N) / (math.factorial(x1) * math.factorial(x2) * math.factorial(x3)) / 3 ** N
    return p, s_obs


def binom_tail(kk, n): return sum(math.comb(n, i) for i in range(kk, n + 1)) / 2 ** n


def analyse(units, label, boot=1000, cert=True):
    R = {'label': label, 'units': len(units)}
    XU, XR = unit_matrices(units)
    uw = XU.sum(0); rw = XR.sum(0)
    w = {s: int(uw[si[s]]) for s in SPL}
    sc = MU @ uw; order = np.argsort(-sc, kind='stable')
    best = UT[order[0]]
    R['unrooted'] = dict(best=sorted(k(s) for s in best), score=int(sc[order[0]]), runner_up=int(sc[order[1]]),
                         runner_up_tree=sorted(k(s) for s in UT[order[1]]),
                         ties_at_top=int((sc == sc[order[0]]).sum()),
                         weights={x: w.get(frozenset(x), 0) for x in ['CG', 'DE', 'CDEFG', 'DEF', 'EF', 'DF', 'CEG', 'BCG', 'CD']})
    false = {s: v for s, v in w.items() if s not in best and v > 0}
    mint = min(w[s] for s in best); fmax = max(false.values())
    R['global_cert'] = dict(min_true=mint, min_true_split=[k(s) for s in best if w[s] == mint],
                            max_false=fmax, max_false_split=[k(s) for s, v in false.items() if v == fmax],
                            margin=mint - fmax, holds=mint > fmax)
    loc = sorted((w[s] - v, k(s), w[s], k(c), v) for c, v in false.items() for s in best if crosses(c, s))
    R['local_cert'] = dict(min_margin=loc[0][0], binding=[dict(true=a, w_true=b, rival=c, w_rival=d) for m, a, b, c, d in loc if m == loc[0][0]],
                           holds=loc[0][0] > 0)
    # per-family rival test (rival_test.py logic: strongest crossing rival, one-sided binomial)
    rv = []
    for s in sorted(best, key=lambda s: -w[s]):
        rivals = sorted(((v, k(c)) for c, v in false.items() if crosses(c, s)), reverse=True)
        v, r = rivals[0]
        rv.append(dict(family=k(s), support=w[s], strongest_rival=r, rival_support=v, p_one_sided=binom_tail(w[s], w[s] + v)))
    R['rival_test'] = rv
    de, ef, df = w[frozenset('DE')], w[frozenset('EF')], w[frozenset('DF')]
    p, s_obs = trinomial_p((de, ef, df))
    R['trio_DE_EF_DF'] = dict(counts=[de, ef, df], p=p, s_obs=s_obs)
    # rooted, min weight 2 (lai.py convention) and min weight 1
    for mn in (2, 1):
        rwm = np.where(rw >= mn, rw, 0); rs = MR @ rwm; o = np.argsort(-rs, kind='stable')
        rwd = {c: int(rwm[ci[c]]) for c in CLU}
        base = [frozenset(x) for x in ['AB', 'CG', 'DE', 'DEF']]
        rt = {}
        for nm, t in [('AB-edge', base + [frozenset('CDEFG')]), ('DEF-edge', base + [frozenset('ABCG')]), ('central', base)]:
            rt[nm] = dict(score=sum(rwd[c] for c in t), clusters={k(c): rwd[c] for c in t})
        # all 11 rootings of the unrooted optimum
        canon = lambda c: c if 'A' not in c else ALL - c
        roots = []
        for a, t in enumerate(RT):
            if frozenset(canon(c) for c in t if len(c) <= 5) == best:
                mx = [c for c in t if not any(c < d for d in t)]
                kids = mx if len(mx) == 2 else mx + [ALL - mx[0]]
                roots.append((int(rs[a]), ' | '.join(sorted((k(c) for c in kids), key=len)), sorted(k(c) for c in t)))
        R[f'rooted_min{mn}'] = dict(best=sorted(k(c) for c in RT[o[0]]), score=int(rs[o[0]]),
                                    runner_up=sorted(k(c) for c in RT[o[1]]), runner_up_score=int(rs[o[1]]),
                                    ties_at_top=int((rs == rs[o[0]]).sum()), named=rt,
                                    rootings_of_unrooted_optimum=sorted(roots, reverse=True))
    # bootstrap stability of the unrooted optimum's splits and of the rooted optimum's clusters
    rng = np.random.default_rng(20261122); n = len(units)
    stab_u = collections.Counter(); stab_r = collections.Counter()
    for _ in range(boot):
        c = np.bincount(rng.integers(0, n, n), minlength=n)
        bu = c @ XU; br = c @ XR; br = np.where(br >= 2, br, 0)
        stab_u.update(UT[int(np.argmax(MU @ bu))]); stab_r.update(RT[int(np.argmax(MR @ br))])
    R['stability_unrooted'] = {k(s): stab_u[s] / boot for s in best}
    R['stability_rooted'] = {k(c): stab_r[c] / boot for c in RT[o[0]]}
    if cert: R['gamma_confidence'] = gamma_cert(XU, best, n)
    return R


def gamma_cert(XU, best, L):
    """Thm 6.9 sample certificate (gamma_confidence.py, Lai part, reused verbatim through exec of its definitions)."""
    src = open(os.path.join(ROOT, 'identifiability', 'code', 'gamma_confidence.py')).read()
    g = {}; exec(src.split("res = {'delta': DELTA}")[0], g)
    X7, S7 = g['all_splits'](7); sol = g['Solver'](X7, S7)
    ix = {c: j for j, c in enumerate(W)}
    tok = lambda s: frozenset(ix[c] for c in s)
    cnt = {tok(s): int(v) for s, v in zip(SPL, XU.sum(0)) if v}
    T = {tok(s) for s in best}
    _, That = sol.best({v: cnt.get(v, 0) / L for v in S7}); assert That == T
    g_lo, g_hat, _ = g['certify'](sol, cnt, L, T, 7)
    scaled = {}
    for kk in (1, 1.5, 2, 3, 4, 5, 6, 8):
        c2 = {v: int(round(kk * x)) for v, x in cnt.items()}
        scaled[kk] = g['certify'](sol, c2, int(round(kk * L)), T, 7)[0]
    first = next((kk for kk, v in scaled.items() if v > 0), None)
    return dict(L=L, gamma_hat=g_hat, gamma_lo=g_lo, scaled=scaled, first_multiple_certified=first)


def majority_reference(units):
    """Sensitivity: re-polarise each unit against its majority reading (a class holding a strict majority of the
    witnesses present; otherwise the printed text is kept) instead of the printed (A-corrected) text. Classes = listed groups + the base (witnesses not listed). F lacks vv. 1-52 (prologue). Other lacunae are not
    modelled, so this is approximate."""
    out = []
    for u in units:
        present = ALL - (frozenset('F') if u['line'] <= 52 else frozenset())
        gs = [frozenset(g) & present for g in u['groups'][1:] if g]
        base = present - frozenset().union(*gs) if gs else present
        classes = [base] + [g for g in gs if g]
        ref = max(classes, key=len)
        if 2 * len(ref) <= len(present) or len(ref) == len(base): ref = base   # no strict majority, or tie: printed text
        dep = [c for c in classes if c is not ref and c]
        out.append(dict(u, groups=[''] + [k(c) for c in dep]))
    return out


if __name__ == '__main__':
    V2, nsplit = build()
    print('units_v2:', len(V2), 'units;', nsplit, 'elementary units from', len(SPLITS), 'composite notes; n056 dropped')
    V1 = U
    RES = {'v1_published': analyse(V1, 'v1 published (814 notes)'),
           'v1_minus_n056': analyse([u for u in V1 if u['page'] not in DROP_PAGES], 'v1 without n056 (no splits)'),
           'v2': analyse(V2, 'v2: n056 dropped, 34 composite notes split')}
    RES['v2_majority_reference'] = analyse(majority_reference(V2), 'v2 against the majority reading', cert=False)
    RES['v2_majority_reference']['note'] = majority_reference.__doc__
    json.dump(RES, open(os.path.join(D, 'results_v2.json'), 'w'), indent=1, ensure_ascii=False, default=str)
    for key, r in RES.items():
        print(f"\n== {r['label']}  units={r['units']}")
        u = r['unrooted']; print(' unrooted', u['best'], u['score'], 'runner-up', u['runner_up'], 'ties', u['ties_at_top'])
        print(' weights', u['weights'])
        print(' global', r['global_cert']); print(' local', r['local_cert'])
        print(' rival', [(x['family'], x['support'], x['strongest_rival'], x['rival_support'], round(x['p_one_sided'], 5)) for x in r['rival_test']])
        print(' trio DE/EF/DF', r['trio_DE_EF_DF'])
        for mn in (2, 1):
            q = r[f'rooted_min{mn}']
            print(f' rooted(min{mn}) best', q['best'], q['score'], 'runner-up', q['runner_up'], q['runner_up_score'], 'ties', q['ties_at_top'])
            print('   named', {a: b['score'] for a, b in q['named'].items()}, q['named']['AB-edge']['clusters'], 'ABCG', q['named']['DEF-edge']['clusters']['ABCG'])
            print('   rootings', [(s, o) for s, o, _ in q['rootings_of_unrooted_optimum']])
        print(' stability unrooted', r['stability_unrooted']); print(' stability rooted', r['stability_rooted'])
        if 'gamma_confidence' in r: print(' gamma', r['gamma_confidence'])
