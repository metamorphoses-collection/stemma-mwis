"""Test 1: is the device's outcome distribution approximately Gibbs over independent sets?
For every existing sampler-test emulation (instance x layout x T x encoding x noise), fit
    P(S) = exp(beta * eps_S) / Z  over a support set Omega,
eps_S = sum_{i in S} delta_i - sum_{i<j in S} V_ij (classical end-of-sweep reward, rad/us, ALL vdW tails, the delta_i the run
used: unweighted delta_f; w<a> delta_f * w_eff(a)). Omega = (i) all independent sets, (ii) maximal sets only.
beta by maximum likelihood of the device distribution restricted to Omega (= min KL(p_dev|Omega || Gibbs)), 1-D.
Noisy runs: fitted on the post-readout exact distribution AND on the readout-mitigated one (known p_fp, p_fn inverted).
Goodness: TV(p_dev|Omega, Gibbs) and, for Omega = maximal, TV over the history shares.
Joint fit: one beta per setting (T, encoding, noise) shared by the three sound-layout instances
(priscus register, lai opt20.0, lawrence opt12.0); TV with the shared beta vs per-instance beta.
Output: calib_20261123/gibbs_fit.json + gibbs_fit.md"""
import json, os, glob
import numpy as np
from scipy.optimize import minimize_scalar
from calib_common import INST, SOUT, OUT, DELTAF, independent_sets, maximal_sets, lab, energies, mitigate, tv
from sampler_test_common import weff


def fit(p, e):
    """p, e arrays over Omega (p renormalised). Returns beta, TV."""
    p = p / p.sum()
    def nll(b):
        s = b * e; s = s - s.max(); q = np.exp(s); q /= q.sum(); return -float((p * np.log(q + 1e-300)).sum())
    r = minimize_scalar(nll, bounds=(-5, 30), method='bounded', options=dict(xatol=1e-5))
    b = float(r.x); s = b * e; s -= s.max(); q = np.exp(s); q /= q.sum()
    return b, tv(p, q), q


def joint_fit(items):
    def nll(b):
        t = 0
        for p, e in items:
            p = p / p.sum(); s = b * e; s = s - s.max(); q = np.exp(s); q /= q.sum(); t -= float((p * np.log(q + 1e-300)).sum())
        return t
    b = float(minimize_scalar(nll, bounds=(-5, 30), method='bounded', options=dict(xatol=1e-5)).x)
    tvs = []
    for p, e in items:
        p = p / p.sum(); s = b * e; s -= s.max(); q = np.exp(s); q /= q.sum(); tvs.append(tv(p, q))
    return b, tvs


res = {}
for f in sorted(glob.glob(os.path.join(SOUT, 'emu_*.json'))):
    E = json.load(open(f)); name = E['instance']; g = INST[name]; F = g['fam']
    P = np.array([E['positions_um'][x] for x in F], float)
    IND = [lab(g, s) for s in independent_sets(g)]; MAX = [lab(g, s) for s in maximal_sets(g)]
    job = os.path.basename(f)[4:-5]; res[job] = {}
    for key, r in E['runs'].items():
        delta = np.full(len(F), DELTAF) if r['a'] is None else np.array([DELTAF * weff(g, r['a'])[x] for x in F])
        eps = energies(name, P, delta)
        dists = {'raw': r['exact']}
        if r['noise'] == 'noisy': dists['mitigated'] = mitigate(r['exact'], name)
        for dk, d in dists.items():
            rec = dict(T=r['T'], mode=r['mode'], noise=r['noise'], dist=dk)
            for om, OM in (('independent', IND), ('maximal', MAX)):
                p = np.array([d.get(h, 0.0) for h in OM]); e = np.array([eps[h] for h in OM])
                b, t, q = fit(p, e)
                rec[om] = dict(beta=b, tv=t, mass=float(p.sum()))
                if om == 'maximal': rec[om]['shares_dev'] = (p / p.sum()).tolist(); rec[om]['shares_gibbs'] = q.tolist(); rec[om]['labels'] = OM
            res[job][f'{key}|{dk}'] = rec
json.dump(res, open(os.path.join(OUT, 'gibbs_fit.json'), 'w'), indent=1)

# ---------------- tables
md = ['# Test 1 — Gibbs fit of the existing sampler-test emulations\n',
      'P(S) ∝ exp(β·ε_S), ε_S = Σ_{i∈S} δ_i − Σ V_ij (rad/µs, all tails). β in µs (1/(rad/µs)). '
      'TV = total-variation distance between device distribution restricted to Ω and the fitted Gibbs law. '
      '"mit" = noisy run after inverting the known readout channel.\n']
SOUND = ['priscus_register', 'lai_opt20.0', 'lawrence_opt12.0']
for job in sorted(res):
    md.append(f'\n## {job}\n')
    md.append('| setting | β (indep.) | TV indep. | β (maximal) | TV maximal | P(maximal) |'); md.append('|---|---|---|---|---|---|')
    for k, r in sorted(res[job].items(), key=lambda t: (t[1]['noise'], t[1]['dist'], t[1]['mode'], t[1]['T'])):
        md.append(f"| {r['mode']} {r['noise']}{' mit' if r['dist'] == 'mitigated' else ''} T{r['T']} | {r['independent']['beta']:.2f} | "
                  f"{r['independent']['tv']:.3f} | {r['maximal']['beta']:.2f} | {r['maximal']['tv']:.3f} | {r['maximal']['mass']:.3f} |")
# joint beta across the three sound-layout instances, per setting
md.append('\n## One β per setting across the three sound-layout instances (priscus register, lai opt20.0, lawrence opt12.0)\n')
md.append('Ω = independent sets. Per-instance β in brackets; TV with shared β (per-instance TV in brackets).\n')
md.append('| setting | shared β | β priscus / lai / lawrence | TV priscus | TV lai | TV lawrence |'); md.append('|---|---|---|---|---|---|')
joint = {}
keys = sorted(set(res['priscus_register']) & set(res['lai_opt20.0']) & set(res['lawrence_opt12.0']),
              key=lambda k: (res['priscus_register'][k]['noise'], res['priscus_register'][k]['dist'], res['priscus_register'][k]['mode'], res['priscus_register'][k]['T']))
for k in keys:
    items = []
    for job in SOUND:
        E = json.load(open(os.path.join(SOUT, f'emu_{job}.json'))); name = E['instance']; g = INST[name]; F = g['fam']
        P = np.array([E['positions_um'][x] for x in F], float); r0 = E['runs'][k.split('|')[0]]
        delta = np.full(len(F), DELTAF) if r0['a'] is None else np.array([DELTAF * weff(g, r0['a'])[x] for x in F])
        eps = energies(name, P, delta); d = r0['exact'] if k.endswith('raw') else mitigate(r0['exact'], name)
        IND = [lab(g, s) for s in independent_sets(g)]
        items.append((np.array([d.get(h, 0.0) for h in IND]), np.array([eps[h] for h in IND])))
    b, tvs = joint_fit(items); r = res['priscus_register'][k]
    per = [res[j][k]['independent'] for j in SOUND]
    joint[k] = dict(beta_shared=b, tv_shared=tvs, beta_each=[x['beta'] for x in per], tv_each=[x['tv'] for x in per])
    md.append(f"| {r['mode']} {r['noise']}{' mit' if r['dist'] == 'mitigated' else ''} T{r['T']} | {b:.2f} | "
              + ' / '.join(f"{x['beta']:.2f}" for x in per) + ' | ' + ' | '.join(f'{t:.3f} ({x["tv"]:.3f})' for t, x in zip(tvs, per)) + ' |')
json.dump(joint, open(os.path.join(OUT, 'gibbs_joint.json'), 'w'), indent=1)
open(os.path.join(OUT, 'gibbs_fit.md'), 'w').write('\n'.join(md) + '\n')
print('\n'.join(md[-40:]))
