"""Summarise bedier_sensitivity.py: per text x scenario, share of replicates whose verdict equals the frozen one; per
scenario, the range of the series counts (supported / contradicted / partly silent) and of the pooled one-sided p of
the editors-vs-random-tops score (bedier_null.py statistic), recomputed on each replicate.
Usage: bedier_sensitivity_summary.py DIR [MAXREP]   (reads DIR/baseline.json, DIR/replicates.jsonl; writes DIR/summary.json, .md)"""
import json, os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bedier_sensitivity import TEXTS, SCEN, SCORE

SHORT = {'isidorus_chronica_maiora': 'Isidore, Chron.', 'claudianus_raptus_proserpinae': 'Claudian, Raptus',
         'prosper_epitoma_chronicon': 'Prosper', 'historia_brittonum': 'Hist. Britt.', 'isidorus_historia_gothorum': 'Isidore, Goth.',
         'fortunatus_vita_martini': 'Vita Martini', 'claudianus_maior': 'Claudian, maiora', 'sidonius_carmina': 'Sidonius carm.',
         'ennodius_opera': 'Ennodius', 'paulus_historia_romana': 'Paul', 'victor_vitensis_historia': 'Victor',
         'jordanes_romana': 'Romana', 'eugippius_vita_severini': 'Eugippius', 'marcellinus_comes_chronicon': 'Marcellinus'}
AB = {'supported': 'S', 'contradicted': 'C', 'partly_silent': 'Si'}


def pooled_p(nulls, obs_scores):
    D = {0: 1.0}
    for nl in nulls:
        nd = collections.defaultdict(float)
        for s0, p0 in D.items():
            for v, q in nl.items(): nd[s0 + SCORE[v]] += p0 * q
        D = nd
    obs = sum(obs_scores); return sum(v for k, v in D.items() if k >= obs - 1e-9)


def main():
    d = sys.argv[1]; maxrep = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
    base = {r['text']: r for r in json.load(open(os.path.join(d, 'baseline.json')))}
    R = collections.defaultdict(dict)
    for line in open(os.path.join(d, 'replicates.jsonl')):
        r = json.loads(line)
        if r['rep'] < maxrep: R[r['scenario']][(r['text'], r['rep'])] = r
    p0 = pooled_p([b['null'] for b in base.values() if b['null']], [SCORE[b['verdict']] for b in base.values() if b['null']])
    out = dict(baseline=dict(counts=collections.Counter(b['verdict'] for b in base.values()), pooled_p=p0), scenarios={})
    for sc in SCEN:
        if sc not in R: continue
        rows = R[sc]; reps = sorted({k[1] for k in rows})
        per = {}
        for t in TEXTS:
            vs = [rows[(t, i)]['verdict'] for i in reps if (t, i) in rows]
            per[t] = dict(n=len(vs), unchanged=sum(v == base[t]['verdict'] for v in vs) / len(vs),
                          dist=dict(collections.Counter(vs)),
                          gate_fail=sum(not rows[(t, i)]['gate_200'] for i in reps if (t, i) in rows) / len(vs))
        cnt, ps, allsame = [], [], 0
        for i in reps:
            if not all((t, i) in rows for t in TEXTS): continue
            vs = [rows[(t, i)]['verdict'] for t in TEXTS]
            c = collections.Counter(vs); cnt.append((c['supported'], c['contradicted'], c['partly_silent']))
            allsame += all(rows[(t, i)]['verdict'] == base[t]['verdict'] for t in TEXTS)
            used = [rows[(t, i)] for t in TEXTS if rows[(t, i)]['null']]
            ps.append(pooled_p([r['null'] for r in used], [SCORE[r['verdict']] for r in used]))
        n = len(cnt); ps.sort()
        rng = lambda k: (min(c[k] for c in cnt), max(c[k] for c in cnt))
        out['scenarios'][sc] = dict(per_text=per, n_complete=n, all_14_unchanged=allsame / n if n else None,
                                    counts_exact_4_4_6=sum(c == (4, 4, 6) for c in cnt) / n if n else None,
                                    supported=rng(0), contradicted=rng(1), silent=rng(2),
                                    supported_mean=sum(c[0] for c in cnt) / n, contradicted_mean=sum(c[1] for c in cnt) / n,
                                    silent_mean=sum(c[2] for c in cnt) / n,
                                    p_median=ps[n // 2], p_95=ps[int(0.95 * (n - 1))], p_max=ps[-1],
                                    share_p_lt_001=sum(p < 1e-3 for p in ps) / n,
                                    non_optimal=sum(r['non_optimal'] for r in rows.values()),
                                    secs=sum(r['secs'] for r in rows.values()))
    json.dump(out, open(os.path.join(d, 'summary.json'), 'w'), indent=1, default=list)
    S = [s for s in SCEN if s in out['scenarios']]
    L = ['| Text (frozen) | ' + ' | '.join(S) + ' |', '|---|' + '---|' * len(S)]
    for t in TEXTS:
        L.append(f"| {SHORT[t]} ({AB[base[t]['verdict']]}) | " + ' | '.join(f"{100 * out['scenarios'][s]['per_text'][t]['unchanged']:.0f}" for s in S) + ' |')
    for lab, k in [('supported (min-max)', 'supported'), ('contradicted', 'contradicted'), ('silent', 'silent')]:
        L.append(f'| {lab} | ' + ' | '.join('%d-%d' % tuple(out['scenarios'][s][k]) for s in S) + ' |')
    L.append('| all 14 unchanged % | ' + ' | '.join(f"{100 * out['scenarios'][s]['all_14_unchanged']:.0f}" for s in S) + ' |')
    L.append('| pooled p, median | ' + ' | '.join(f"{out['scenarios'][s]['p_median']:.1e}" for s in S) + ' |')
    L.append('| pooled p, max | ' + ' | '.join(f"{out['scenarios'][s]['p_max']:.1e}" for s in S) + ' |')
    txt = f"Baseline: {dict(out['baseline']['counts'])}, pooled p = {p0:.2e}\n\n" + '\n'.join(L) + '\n'
    open(os.path.join(d, 'summary.md'), 'w').write(txt); print(txt)
    for s in S:
        for t in TEXTS:
            pt = out['scenarios'][s]['per_text'][t]
            if pt['unchanged'] < 0.95 or pt['gate_fail'] > 0: print(s, t, pt)


if __name__ == '__main__':
    main()
