"""EXPLORATORY (not pre-registered; written after the classification of 2026-11-13): are editors' top divisions
better supported than random divisions of the same witnesses?
For each classified text, the stable families of its frozen classification run (classification_a10/<text>.json) are
fixed; they do not depend on the editor's top. Every pair of disjoint branches drawn from the analysis witnesses with the
editor's branch sizes (witnesses outside his top stay outside) is scored by the frozen rule (bedier_top.py): a branch with >= 2 witnesses is supported if it is a
stable family, contradicted if a stable family of weight >= 5 crosses it; the split is supported / contradicted /
partly silent as in the protocol. Splits are enumerated exactly.
Per text: rank of the editor's split among all size-matched splits (score: supported 2, partly silent 1, contradicted 0),
one-sided p = P(random score >= editor score). Pooled: sum of editor scores against the exact null distribution of the
sum (independent texts), and the same for 'supported' and 'not contradicted' counts.
Usage: bedier_null.py text_id ... --out FILE"""
import json, itertools, sys, os, argparse, collections
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'historic', 'bedier')
CROSS_W = 5


def crosses(a, b): return bool(a & b) and not (a <= b or b <= a)


def verdict(branches, stable):
    """stable: list of (frozenset, weight)."""
    st = []
    S = {f for f, _ in stable}
    for b in branches:
        if len(b) < 2: continue
        if b in S: st.append('supported')
        elif any(w >= CROSS_W and crosses(f, b) for f, w in stable): st.append('contradicted')
        else: st.append('silent')
    if not st: return None
    return 'contradicted' if 'contradicted' in st else ('supported' if all(s == 'supported' for s in st) else 'partly_silent')


SCORE = {'supported': 2, 'partly_silent': 1, 'contradicted': 0}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('texts', nargs='+'); ap.add_argument('--out', required=True); a = ap.parse_args()
    res = {}
    for t in a.texts:
        d = json.load(open(os.path.join(B, 'classification_a10', t + '.json')))
        W = list(d['analysis_witnesses'])
        top = [frozenset(x) for x in d['editor_top']]
        stable = [(frozenset(k), w) for k, w, s in d['stable_families']]
        ed = verdict(top, stable)
        assert ed == d['classification'], (t, ed, d['classification'])
        if len(top) != 2:
            res[t] = dict(skipped='editor top has one branch among the analysis witnesses', editor=ed); continue
        # random disjoint branches of the editor's sizes (witnesses outside his top stay outside, as in his top)
        s1, s2 = sorted([len(top[0]), len(top[1])]); n = len(W)
        dist = collections.Counter(); seen = set()
        for c1 in itertools.combinations(W, s1):
            rest = [w for w in W if w not in c1]
            for c2 in itertools.combinations(rest, s2):
                key = frozenset([frozenset(c1), frozenset(c2)])
                if key in seen: continue
                seen.add(key); dist[verdict([frozenset(c1), frozenset(c2)], stable)] += 1
        N = sum(dist.values()); e = SCORE[ed]
        p = sum(v for v_, v in dist.items() if SCORE[v_] >= e) / N
        res[t] = dict(n_witnesses=n, branch_sizes=[s1, s2], editor=ed, n_splits=N,
                      null={kk: v / N for kk, v in dist.items()}, p_one_sided=p)
    # pooled exact null of the score sum, and of counts
    used = [r for r in res.values() if 'null' in r]
    def conv(fun):
        D = {0: 1.0}
        for r in used:
            nd = collections.defaultdict(float)
            for s0, p0 in D.items():
                for v, q in r['null'].items(): nd[s0 + fun(v)] += p0 * q
            D = nd
        return D
    out = {}
    for name, fun in [('score_sum', lambda v: SCORE[v]), ('n_supported', lambda v: v == 'supported'), ('n_not_contradicted', lambda v: v != 'contradicted')]:
        D = conv(fun); obs = sum(fun(r['editor']) for r in used)
        out[name] = dict(observed=obs, null_mean=sum(k * v for k, v in D.items()), p_one_sided=sum(v for k, v in D.items() if k >= obs))
    json.dump(dict(texts=res, pooled=out, n_texts=len(used)), open(a.out, 'w'), indent=1)
    for t, r in res.items():
        if 'null' in r: print(f"{t:30s} n={r['n_witnesses']:2d} sizes={r['branch_sizes']} splits={r['n_splits']:6d} editor={r['editor']:13s} "
                              f"null sup={r['null'].get('supported',0):.3f} sil={r['null'].get('partly_silent',0):.3f} con={r['null'].get('contradicted',0):.3f} p={r['p_one_sided']:.3f}")
        else: print(f"{t:30s} {r['skipped']} ({r['editor']})")
    print('POOLED', json.dumps(out))


if __name__ == '__main__':
    main()
