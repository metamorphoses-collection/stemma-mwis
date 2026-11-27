"""CHECK (2026-11-23): collects the review outputs into one table (data/bedier/checks_20261123/summary.json):
per-text verdicts under every condition, B25/B29 family counts on the frozen run, and the pooled statistic for the B10
compatible-division null (exact convolution over texts with a non-empty compatible set; one-sided P(sum >= observed))."""
import json, os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bedier_sensitivity import TEXTS, ROOT, SCORE, crosses
from rooted_mwis import weigh
O = os.path.join(ROOT, 'checks_20261123')
C = ['printed', 'maj_all', 'maj_span', 'nosp_strict', 'nosp_broad', 'nosp_union']
R = {c: {r['text']: r for r in json.load(open(os.path.join(O, f'cond_{c}.json')))['texts']} for c in C}
P = {c: json.load(open(os.path.join(O, f'cond_{c}.json')))['pooled'] for c in C}
comp = json.load(open(os.path.join(O, 'compat_null.json')))
rows = []
for t in TEXTS:
    d = json.load(open(os.path.join(ROOT, 'classification_a10', t + '.json')))
    top = [frozenset(b) for b in d['editor_top']]; br = [b for b in top if len(b) >= 2]
    st = [(frozenset(k), w) for k, w, s in d['stable_families']]
    U = json.load(open(os.path.join(ROOT, 'runs', t, 'units.json'))); Ws = frozenset(d['analysis_witnesses'])
    F = []
    for u in U:
        fs = [frozenset(x for x in g if x in Ws) for g in u['groups'][1:]]
        fs = [g for g in fs if len(g) >= 2 and g != Ws]
        if fs: F.append(fs)
    cand = weigh(F, 2)
    rows.append(dict(text=t, editor_top=d['editor_top'], n_analysis_witnesses=len(d['analysis_witnesses']),
                     candidate_families_weight_ge2=len(cand), families_in_optimum=len(d['mwis']), stable_families=len(st),
                     stable_weight_ge5=sum(1 for f, w in st if w >= 5),
                     stable_ge5_crossing_a_classified_branch=sum(1 for f, w in st if w >= 5 and any(crosses(f, b) for b in br)),
                     crossing=[''.join(sorted(f)) + f' {w}' for f, w in st if w >= 5 and any(crosses(f, b) for b in br)],
                     classified_branches=len(br),
                     verdict={c: R[c][t]['verdict'] for c in C},
                     majority_branch_over_half=[''.join(sorted(b)) for b in br if 2 * len(b) > len(d['analysis_witnesses'])]))
# B10 pooled: compatible ('all' stable families) null
pool = {}
for which in ('all', 'nonmax'):
    D = {0: 1.0}; obs = 0; used = []; dropped = []
    for t in TEXTS:
        r = comp[t][which]
        if not r: continue
        if not r['n_compatible']: dropped.append(t); continue
        used.append(t); obs += SCORE[r['editor']]
        nd = collections.defaultdict(float)
        for s0, p0 in D.items():
            for v, q in r['shares'].items(): nd[s0 + SCORE[v]] += p0 * q
        D = nd
    pool[which] = dict(texts=len(used), dropped_no_compatible_division=dropped, observed_score_sum=obs,
                       null_mean=sum(k * v for k, v in D.items()), null_min=min(D), p_one_sided=sum(v for k, v in D.items() if k >= obs),
                       editor_not_compatible=[t for t in used if not comp[t][which]['editor_compatible']],
                       median_share_contradicted=sorted(comp[t][which]['shares'].get('contradicted', 0) for t in used)[len(used) // 2])
out = dict(per_text=rows, pooled_random_null={c: P[c] for c in C}, pooled_compatible_null=pool)
json.dump(out, open(os.path.join(O, 'summary.json'), 'w'), indent=1)
for r in rows: print(r['text'][:28].ljust(28), r['n_analysis_witnesses'], r['candidate_families_weight_ge2'], r['families_in_optimum'], r['stable_families'], r['stable_weight_ge5'], r['stable_ge5_crossing_a_classified_branch'], r['crossing'], r['majority_branch_over_half'])
print(json.dumps(pool, indent=1))
