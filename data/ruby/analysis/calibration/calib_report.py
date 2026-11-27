"""Test 3c/4: compare device shares with bootstrap shares across all graded cases.
Inputs: calib_20261123/cases.json, calibration.json, runs.json. Outputs: calib_20261123/report_tables.md, scatter.tsv,
summary.json, scatter_T6000.png."""
import json, os
import numpy as np
from calib_common import OUT, tv

C = {c['case']: c for c in json.load(open(os.path.join(OUT, 'cases.json')))}
R = json.load(open(os.path.join(OUT, 'runs.json')))
CAL = json.load(open(os.path.join(OUT, 'calibration.json')))
ENC = ['raw', 'th', 'th0', 'lai0', 'kl0']
SETS = [(T, n, r) for T in (6000, 4000, 2000) for n, r in (('ideal', 'ideal'), ('noisy', 'raw'), ('noisy', 'mit'))]
md, summ, rows = [], {}, ['case\thistory\tbootstrap\tsetting\tencoding\tdevice_exact\tdevice_samples']
OUTOF = {'lai0': {'lai_full'}}   # in-sample: beta_lai was fitted on lai_full


def metrics(cases, T, n, r, e, samples=False):
    b, d, tvs, top_err, pm = [], [], [], [], []
    for cn in cases:
        k = f'T{T}|{n}|{r}|{e}'
        if k not in R.get(cn, {}): continue
        x = R[cn][k]; c = C[cn]
        s = np.array(x['share_used'] if not samples else (x['share_samples_mitigated'] if r == 'mit' else x['share_samples']))
        p = np.array(c['bootstrap']); b += p.tolist(); d += s.tolist(); tvs.append(tv(s, p))
        i = int(np.argmax(p)); top_err.append(abs(s[i] - p[i])); pm.append(x['P_maximal'])
    if not tvs: return None
    b, d = np.array(b), np.array(d)
    return dict(n_cases=len(tvs), mean_tv=float(np.mean(tvs)), max_tv=float(np.max(tvs)), mae_all=float(np.abs(b - d).mean()),
                mae_top=float(np.mean(top_err)), r=float(np.corrcoef(b, d)[0, 1]), mean_P_maximal=float(np.mean(pm)))


allc = list(C)
md.append('# Calibrated encoding vs bootstrap: summary over all 16 cases (exact device distributions)\n')
md.append('Encodings: raw = preps\' w_eff(0.3) control; th = scale from thermometer only; th0 = thermometer + register zero point; '
          'lai0 = β from Lai-full curve + zero point; kl0 = ORACLE scale θ_KL (fitted to the bootstrap) + thermometer β + zero point.\n')
md.append('| setting | encoding | mean TV | max TV | MAE top share | MAE all shares | Pearson r | mean P(maximal) |')
md.append('|---|---|---|---|---|---|---|---|')
for T, n, r in SETS:
    for e in ENC:
        m = metrics(allc, T, n, r, e)
        if m is None: continue
        summ[f'T{T}|{n}|{r}|{e}'] = m
        md.append(f"| T{T} {n}{'/'+r if n == 'noisy' else ''} | {e} | {m['mean_tv']:.3f} | {m['max_tv']:.3f} | {m['mae_top']:.3f} | "
                  f"{m['mae_all']:.3f} | {m['r']:.3f} | {m['mean_P_maximal']:.3f} |")
# encoding ceiling
ce = dict(pred=np.mean([c['tv_softmax_pred'] for c in C.values()]), kl=np.mean([c['tv_softmax_kl'] for c in C.values()]),
          mvn=np.mean([c['tv_mvn'] for c in C.values()]))
md.append(f"\nEncoding ceiling (a perfect Gibbs device, no physics): mean TV(softmax(θ_pred W), bootstrap) = {ce['pred']:.3f}; "
          f"with θ_KL = {ce['kl']:.3f}; Gaussian (multinomial-probit) approximation = {ce['mvn']:.3f}.\n")
summ['ceiling'] = ce

# 4000-sample version for the main encodings
md.append('\n## Same with 4000-sample shares (finite-shot estimate; mitigated reads invert the readout channel on the counts)\n')
md.append('| setting | encoding | mean TV | MAE top | Pearson r |'); md.append('|---|---|---|---|---|')
for T, n, r in SETS:
    for e in ('raw', 'th0', 'lai0'):
        m = metrics(allc, T, n, r, e, samples=True)
        if m: md.append(f"| T{T} {n}{'/'+r if n == 'noisy' else ''} | {e} | {m['mean_tv']:.3f} | {m['mae_top']:.3f} | {m['r']:.3f} |")
        if m: summ[f'samples|T{T}|{n}|{r}|{e}'] = m

# by instance family
md.append('\n## By instance (T6000; mean TV)\n')
md.append('| setting | encoding | Lawrence (7) | Lai (8) | Priscus (1) |'); md.append('|---|---|---|---|---|')
for n, r in (('ideal', 'ideal'), ('noisy', 'raw'), ('noisy', 'mit')):
    for e in ENC:
        cells = []
        for fam in ('lawrence', 'lai', 'priscus'):
            m = metrics([k for k in allc if C[k]['instance'] == fam], 6000, n, r, e); cells.append(f"{m['mean_tv']:.3f}" if m else '')
        md.append(f"| T6000 {n}{'/'+r if n == 'noisy' else ''} | {e} | " + ' | '.join(cells) + ' |')

# per-case table, T6000
for n, r in (('ideal', 'ideal'), ('noisy', 'mit'), ('noisy', 'raw')):
    md.append(f"\n## Per case, T6000 {n}{'/'+r if n == 'noisy' else ''}: top-history share (bootstrap → device), TV in brackets\n")
    md.append('| case | units | z | bootstrap top | raw | th | th0 | lai0 | kl0 | th0: θ_eff/θ | th0 ground = top (gap) | P(max) th0 |')
    md.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
    for cn, c in C.items():
        i = int(np.argmax(c['bootstrap'])); cells = []
        for e in ENC:
            x = R[cn].get(f'T6000|{n}|{r}|{e}')
            cells.append('' if x is None else f"{x['share_used'][i]:.2f} ({tv(x['share_used'], c['bootstrap']):.2f})")
        x = R[cn].get(f'T6000|{n}|{r}|th0'); m = x['meta']
        md.append(f"| {cn} | {c['units']} | {c['z_top_runner']:.2f} | {c['bootstrap'][i]:.3f} | " + ' | '.join(cells) +
                  f" | {m['theta_eff'] / m['theta']:.2f}{' sat' if x['saturated'] else ''} | {'yes' if x['ground']['is_top'] else 'NO ' + x['ground']['state']} ({x['ground']['gap']:.2f}) | {x['P_maximal']:.2f} |")

# scatter
for cn, c in C.items():
    for key, x in R[cn].items():
        for h, pb, pd, ps in zip(c['labels'], c['bootstrap'], x['share_used'],
                                 x.get('share_samples_mitigated') if x['read'] == 'mit' else x['share_samples']):
            rows.append(f"{cn}\t{h}\t{pb:.4f}\t{key.rsplit('|', 1)[0]}\t{key.rsplit('|', 1)[1]}\t{pd:.4f}\t{ps:.4f}")
open(os.path.join(OUT, 'scatter.tsv'), 'w').write('\n'.join(rows) + '\n')
json.dump(summ, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)
open(os.path.join(OUT, 'report_tables.md'), 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))

try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.4), sharex=True, sharey=True)
    for a, (n, r, tt) in zip(ax, (('ideal', 'ideal', 'ideal'), ('noisy', 'raw', 'noisy, raw read'), ('noisy', 'mit', 'noisy, readout-mitigated'))):
        for e, col, mk in (('raw', '#999999', 'o'), ('th0', '#1f5fa8', 's'), ('lai0', '#d9822b', '^')):
            xs, ys = [], []
            for cn, c in C.items():
                x = R[cn].get(f'T6000|{n}|{r}|{e}')
                if x: xs += c['bootstrap']; ys += x['share_used']
            a.scatter(xs, ys, s=16, c=col, marker=mk, label=e, alpha=0.8)
        a.plot([0, 1], [0, 1], 'k:', lw=1); a.set_title(f'T = 6000 ns, {tt}'); a.set_xlabel('bootstrap share')
    ax[0].set_ylabel('device share among maximal outcomes'); ax[0].legend(frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'scatter_T6000.png'), dpi=130)
except Exception as ex:
    print('plot skipped', ex)
