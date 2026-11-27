"""Test 3a: measure beta_eff per device setting, three ways.
Settings: T in {2000, 4000, 6000} ns (on-device), noise in {ideal, noisy}; noisy read both raw (post-readout) and
readout-mitigated (known p_fp/p_fn inverted). Drive as the preps (sampler_test_emulate.py); per-atom detunings from (tail-compensated)
calib_encode.encode_qp (true history rewards = s*W_h + c, all tails; detunings closest to delta_f).
 (1) calibration instance Lawrence full (layout opt12.0): grid of scales s, nominal Delta_eps(E,R) = 3 s in
     {0, .25, .5, 1, 1.5, 2, 2.5, 3.5} rad/us; device shares among maximal outcomes vs the TRUE end-of-sweep reward
     difference (all tails). beta_L = ML Gibbs fit over the maximal histories pooled over the grid (no intercept);
     also logistic fit with intercept (bias b0) and the full curve logit(Delta_eps) for inversion.
 (2) second calibration instance Lai full (opt20.0): grid Delta_eps(top, runner-up) = 18 s in {0, .25, .5, 1, 1.5, 2, 2.5};
     beta_lai = pooled ML Gibbs fit over the 4 maximal histories. (Stability of beta across instances.)
 (3) thermometer (Theo, point 4): an isolated blockaded PAIR (5.5 um apart), detunings 0.65 delta_f +/- Delta_ref/2,
     Delta_ref in {0.5, 1, 2, 3}; beta_th = logit P(a | exactly one excited) / Delta_ref. In emulation a far copy does not
     interact and the noise is per atom, so running the pair separately at the same setting is exactly equivalent to
     carrying it in the same register.
Exact distributions + 4000 samples kept. Output: calib_20261123/calibration.json (rewritten after each block)."""
import json, os, time, sys, zlib
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from calib_common import (INST, OUT, DELTAF, UNITS, X_of, M_of, layout, energies, emulate, mitigate, shares, maximal_sets, lab)
from calib_encode import encode_qp

TS = [2000, 4000, 6000]
INST['thermo'] = dict(fam=['a', 'b'], W={'a': 1, 'b': 1}, conflicts=[('a', 'b')])
P_TH = np.array([[-2.75, 0.0], [2.75, 0.0]])
outp = os.path.join(OUT, 'calibration.json')
res = json.load(open(outp)) if os.path.exists(outp) else {'runs': {}}


def record(key, name, P, what, T, noise, extra):
    if key in res['runs']: return res['runs'][key]
    t0 = time.time(); r = emulate(name, P, what, T, noise == 'noisy', seed=zlib.crc32(key.encode()) % 2**31)
    ML = [lab(INST[name], s) for s in maximal_sets(INST[name])]
    eps = energies(name, P, DELTAF * np.asarray(what))
    rec = dict(instance=name, T=T, noise=noise, what=list(map(float, what)), labels=ML, eps={h: eps[h] for h in ML}, **extra,
               exact=r['exact'], counts=r['counts'])
    rec['share_exact'], rec['P_maximal'] = [x.tolist() if hasattr(x, 'tolist') else x for x in shares(r['exact'], ML)]
    cm = sum(r['counts'].get(h, 0) for h in ML); rec['share_samples'] = [r['counts'].get(h, 0) / cm for h in ML]
    if noise == 'noisy':
        rec['share_mitigated'] = shares(mitigate(r['exact'], name), ML)[0].tolist()
        tot = sum(r['counts'].values()); rec['share_samples_mitigated'] = shares(mitigate({h: v / tot for h, v in r['counts'].items()}, name), ML)[0].tolist()
    res['runs'][key] = rec; json.dump(res, open(outp, 'w'), indent=1)
    print(key, f'{time.time() - t0:.1f}s', [round(x, 3) for x in rec['share_exact']], 'Pmax', round(rec['P_maximal'], 3), flush=True)
    return rec


for name, grid, Dtop in (('lawrence', [0, .25, .5, 1, 1.5, 2, 2.5, 3.5], 3), ('lai', [0, .25, .5, 1, 1.5, 2, 2.5], 18), ('priscus', [0, .25, .5, 1, 1.5, 2, 2.5, 3.5], 41)):
    U = UNITS[name](); w = dict(zip(INST[name]['fam'], X_of(name, U).sum(0).tolist())); P = layout(name)
    for T in TS:
        for noise in ('ideal', 'noisy'):
            for de in grid:
                what, s_used, info = encode_qp(name, w, de / Dtop, P)
                record(f'{name}|T{T}|{noise}|de{de}', name, P, what, T, noise,
                       dict(kind='cal', de_nominal=de, s_req=de / Dtop, s_used=s_used, saturated=info['saturated'], delta=info['delta']))
for T in TS:
    for noise in ('ideal', 'noisy'):
        for dr in (0.5, 1.0, 2.0, 3.0):
            c = 0.65 * DELTAF
            record(f'thermo|T{T}|{noise}|dr{dr}', 'thermo', P_TH, np.array([c + dr / 2, c - dr / 2]) / DELTAF, T, noise,
                   dict(kind='thermo', d_ref=dr))


# ---------------------------------------------------------------- fits
def gibbs_nll(b, items, b0=0.0):
    t = 0.0
    for p, e in items:
        s = b * e; s = s - s.max(); q = np.exp(s); q /= q.sum(); t -= float((np.asarray(p) * np.log(q + 1e-300)).sum())
    return t


def share_key(read):
    return {'ideal': 'share_exact', 'raw': 'share_exact', 'mit': 'share_mitigated'}[read]


fits = {}
for T in TS:
    for noise, read in (('ideal', 'ideal'), ('noisy', 'raw'), ('noisy', 'mit')):
        st = f'T{T}|{noise}|{read}'; fits[st] = {}
        for name in ('lawrence', 'lai', 'priscus'):
            R = [r for k, r in res['runs'].items() if r.get('kind') == 'cal' and r['instance'] == name and r['T'] == T and r['noise'] == noise]
            R.sort(key=lambda r: r['de_nominal'])
            items = [(np.array(r[share_key(read)]), np.array([r['eps'][h] for h in r['labels']])) for r in R]
            b = float(minimize_scalar(lambda x: gibbs_nll(x, items), bounds=(-5, 30), method='bounded').x)
            curve = []
            for r, (p, e) in zip(R, items):
                o = np.argsort(-e); top, nxt = o[0], o[1]
                curve.append(dict(de_nominal=r['de_nominal'], de_true=float(e[top] - e[nxt]), s_used=r['s_used'], saturated=r['saturated'],
                                  logit_top_vs_next=float(np.log(p[top] / p[nxt])) if p[top] > 0 and p[nxt] > 0 else None,
                                  share_top=float(p[top]), P_maximal=r['P_maximal'], shares=p.tolist()))
            rec = dict(beta=b, curve=curve)
            if name == 'lawrence':   # logistic with intercept: logit = b0 + beta * de
                x = np.array([c['de_true'] for c in curve]); y = np.array([c['logit_top_vs_next'] for c in curve], float)
                A = np.vstack([np.ones_like(x), x]).T; coef = np.linalg.lstsq(A, y, rcond=None)[0]
                rec['intercept_fit'] = dict(b0=float(coef[0]), beta=float(coef[1]))
            fits[st][name] = rec
        th = []
        for r in sorted([r for r in res['runs'].values() if r.get('kind') == 'thermo' and r['T'] == T and r['noise'] == noise], key=lambda r: r['d_ref']):
            p = r[share_key(read)]; lg = float(np.log(p[0] / p[1])); th.append(dict(d_ref=r['d_ref'], p_a=p[0], logit=lg, beta=lg / r['d_ref'],
                                                                                   P_one=r['P_maximal']))
        fits[st]['thermo'] = th
res['fits'] = fits
json.dump(res, open(outp, 'w'), indent=1)
for st, f in fits.items():
    print(st, 'beta_lawrence', round(f['lawrence']['beta'], 3), 'intercept-fit', {k: round(v, 3) for k, v in f['lawrence']['intercept_fit'].items()},
          'beta_lai', round(f['lai']['beta'], 3), 'thermo', [(t['d_ref'], round(t['beta'], 3)) for t in f['thermo']])
    print('   lawrence curve', [(round(c['de_true'], 2), None if c['logit_top_vs_next'] is None else round(c['logit_top_vs_next'], 2)) for c in f['lawrence']['curve']])
    print('   lai curve     ', [(round(c['de_true'], 2), None if c['logit_top_vs_next'] is None else round(c['logit_top_vs_next'], 2)) for c in f['lai']['curve']])
# (first pass with LP-vertex detunings kept as calib_20261123/calibration_LPvertex.json: its s = 0 Lawrence point sat on a
#  different LP vertex and flipped sign; this pass uses the continuous least-squares encoding encode_qp.)
