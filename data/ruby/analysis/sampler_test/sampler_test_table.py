"""Task 4: compare emulated shares of the competing maximal histories with the evidence (bootstrap) and with set size.
Share of history h = P(outcome == h) / P(outcome is one of the maximal sets); exact (from the final state + readout
channel) with the 4000-sample estimate alongside. Also P(max-weight set), P(any maximal), P(valid).
Reference models for 'what size alone gives': uniform over maximal sets; fugacity lambda^|S| (fit to the unweighted share).
Output: sampler_test_20261123/tables.md + tables.json"""
import json, os, math, itertools as it
from sampler_test_common import INST, OUT, maximal_sets, independent_sets, lab, conf

B = json.load(open(os.path.join(OUT, 'bootstrap.json')))
BOOTKEY = {'priscus': 'priscus', 'lai': 'lai_v1_register', 'lawrence': 'lawrence'}
md, js = [], {}
JOBS = [(nm, lay, suf) for nm, lay in (('priscus', 'register'), ('lai', 'register'), ('lai', 'opt20.0'), ('lawrence', 'register'),
                                        ('lawrence', 'opt12.0')) for suf in ('', '_longT')]
for name, layout, suf in JOBS:
    p = os.path.join(OUT, f'emu_{name}_{layout}{suf}.json')
    if not os.path.exists(p): continue
    E = json.load(open(p)); g = INST[name]; mx = maximal_sets(g); ML = [lab(g, s) for s in mx]
    ind = {lab(g, s) for s in independent_sets(g)}
    W = {lab(g, s): sum(g['W'][x] for x in s) for s in mx}; SZ = {lab(g, s): len(s) for s in mx}
    order = sorted(ML, key=lambda h: -W[h]); top = order[0]
    bs = B[BOOTKEY[name]]['maximal']
    runs = E['runs']; Ts = sorted({r['T'] for r in runs.values()})
    tab = {}
    for key, r in runs.items():
        ex, ct = r['exact'], r['counts']
        pm = sum(ex.get(h, 0) for h in ML); cm = sum(ct.get(h, 0) for h in ML)
        tab[key] = dict(T=r['T'], mode=r['mode'], noise=r['noise'],
                        share_exact={h: ex.get(h, 0) / pm for h in ML}, share_samples={h: ct.get(h, 0) / cm if cm else float('nan') for h in ML},
                        n_maximal_samples=cm, P_any_maximal=pm, P_top=ex.get(top, 0),
                        P_valid=sum(v for h, v in ex.items() if h in ind), tv=r['tv_samples_vs_exact'],
                        top_outcomes=sorted(ex.items(), key=lambda t: -t[1])[:4])
    js[f'{name}_{layout}{suf}'] = dict(order=order, weight=W, size=SZ, bootstrap={h: bs[h]['bootstrap_share'] for h in ML},
                                  ground_checks=E['ground_checks'], geometry=E['geometry'], settings=tab)
    md.append(f"\n## {name} ({layout} layout{', OFF-DEVICE durations' if suf else ''})\n")
    md.append('ground state (Omega=0, all tails): ' + '; '.join(f"{a}: {v['ground']} {'= MWIS' if v['is_mwis'] else 'NOT MWIS'} gap {v['gap']:.2f}"
                                                                for a, v in E['ground_checks'].items()))
    for noise in ('ideal', 'noisy'):
        md.append(f"\n### {noise}: share among maximal outcomes (exact; 4000-sample value in brackets for the longest T)\n")
        cols = [(m, T) for m in ('unweighted', 'w0.3', 'w0.1', 'w0.0') for T in Ts]
        md.append('| history | size | weight | bootstrap | ' + ' | '.join(f'{m} T{T}' for m, T in cols) + ' |')
        md.append('|' + '---|' * (4 + len(cols)))
        for h in order:
            cells = []
            for m, T in cols:
                t = tab.get(f'T{T}_{m}_{noise}')
                cells.append('' if t is None else f"{t['share_exact'][h]:.3f}" + (f" ({t['share_samples'][h]:.3f})" if T == max(Ts) else ''))
            md.append(f"| {h} | {SZ[h]} | {W[h]} | {bs[h]['bootstrap_share']:.3f} | " + ' | '.join(cells) + ' |')
        cells = []
        for m, T in cols:
            t = tab.get(f'T{T}_{m}_{noise}')
            cells.append('' if t is None else f"{t['P_any_maximal']:.3f}")
        md.append('| P(any maximal) | | | | ' + ' | '.join(cells) + ' |')
        cells = []
        for m, T in cols:
            t = tab.get(f'T{T}_{m}_{noise}')
            cells.append('' if t is None else f"{t['P_top']:.3f}")
        md.append('| P(max-weight set) | | | | ' + ' | '.join(cells) + ' |')
    if name == 'lawrence':
        md.append('\nSize-only reference for E (4 families) vs R (3): uniform 0.5; greedy random order 2/3; fugacity lambda/(1+lambda).')
open(os.path.join(OUT, 'tables.md'), 'w').write('\n'.join(md) + '\n')
json.dump(js, open(os.path.join(OUT, 'tables.json'), 'w'), indent=1, default=float)
print('\n'.join(md))


# ---- compact summary: share of the max-weight history, and whether the share ORDER follows the weight order
def spearman(x, y):
    import numpy as np
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1]) if len(x) > 2 else float('nan')


sm = ['\n## Summary: share of the max-weight history among maximal outcomes (exact) · [Spearman(weight, share) over all maximal histories]\n']
for job, J in js.items():
    Ts = sorted({t['T'] for t in J['settings'].values()}); top = J['order'][0]
    sm.append(f"\n{job}: max-weight = {top}; bootstrap share {J['bootstrap'][top]:.3f}")
    sm.append('| mode | noise | ' + ' | '.join(f'T{T}' for T in Ts) + ' |'); sm.append('|' + '---|' * (2 + len(Ts)))
    for mode in ('unweighted', 'w0.3', 'w0.1', 'w0.0'):
        for noise in ('ideal', 'noisy'):
            cells = []
            for T in Ts:
                t = J['settings'].get(f'T{T}_{mode}_{noise}')
                if t is None: cells.append(''); continue
                rho = spearman([J['weight'][h] for h in J['order']], [t['share_exact'][h] for h in J['order']])
                cells.append(f"{t['share_exact'][top]:.3f}" + ('' if math.isnan(rho) else f" [{rho:+.1f}]"))
            sm.append(f'| {mode} | {noise} | ' + ' | '.join(cells) + ' |')
open(os.path.join(OUT, 'tables.md'), 'a').write('\n'.join(sm) + '\n')
print('\n'.join(sm))
