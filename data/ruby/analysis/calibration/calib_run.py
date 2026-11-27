"""Test 3b: encode every graded case with the calibrated scaling and emulate; control = raw-weight encoding.
Per device setting (T in {2000, 4000, 6000}; ideal | noisy read raw | noisy read readout-mitigated) the calibration
constants come from calib_20261123/calibration.json (measured BEFORE, on other runs):
  beta_th  = thermometer pair, Delta_ref = 1 rad/us (Theo point 4);
  beta_lai = pooled Gibbs fit on the Lai-full calibration curve (calibrate on one instance, apply to all);
  b_h      = ZERO POINT of the register: log device share of history h with all history rewards equal (s = 0 run of
             the same register and setting; weight-independent, so a register property, not a data property).
Encodings (detunings from calib_encode.encode_qp unless 'raw'):
  raw    control: the preps' w_eff = 0.3 + 0.7 w/w_max (case weights), no tail compensation, no calibration;
  th     s = theta_pred / beta_th, theta_pred = 1.702 / sigma_D(top, runner-up)  -- scale only;
  th0    same + zero point;
  lai0   s = theta_pred / beta_lai + zero point (beta from another instance's curve);
  kl0    s = theta_KL / beta_th + zero point; theta_KL fitted to the bootstrap itself (ORACLE scale: isolates device error
         from the softmax-vs-probit encoding error).
For a noisy setting read raw, raw-read constants are used; for mitigated, mitigated-read constants (separate runs).
Exact distribution + 4000 samples per run kept. Output: calib_20261123/runs.json (rewritten after every case)."""
import json, os, time, zlib, sys
import numpy as np
from calib_common import INST, OUT, DELTAF, layout, energies, emulate, mitigate, shares, maximal_sets, lab, tv
from calib_encode import encode_qp, raw_encoding

TS = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '6000,4000,2000').split(',')]
CASES = json.load(open(os.path.join(OUT, 'cases.json')))
CAL = json.load(open(os.path.join(OUT, 'calibration.json')))
outp = os.path.join(OUT, 'runs.json')
res = json.load(open(outp)) if os.path.exists(outp) else {}
READS = (('ideal', 'ideal'), ('noisy', 'raw'), ('noisy', 'mit'))
SK = {'ideal': 'share_exact', 'raw': 'share_exact', 'mit': 'share_mitigated'}


def constants(T, noise, read):
    f = CAL['fits'][f'T{T}|{noise}|{read}']
    bth = [t for t in f['thermo'] if t['d_ref'] == 1.0][0]['beta']
    zero = {}
    for name in ('lawrence', 'lai', 'priscus'):
        r = CAL['runs'][f'{name}|T{T}|{noise}|de0']; p = np.array(r[SK[read]]); b = np.log(np.clip(p, 1e-6, None)); b -= b.mean()
        zero[name] = dict(zip(r['labels'], b.tolist()))
    return dict(beta_th=bth, beta_lai=f['lai']['beta'], zero=zero)


def run(key, name, P, what, T, noise):
    r = emulate(name, P, what, T, noise == 'noisy', seed=zlib.crc32(key.encode()) % 2**31)
    ML = [lab(INST[name], s) for s in maximal_sets(INST[name])]
    out = dict(exact=r['exact'], counts=r['counts'])
    out['share_exact'], out['P_maximal'] = shares(r['exact'], ML); out['share_exact'] = out['share_exact'].tolist()
    cm = sum(r['counts'].get(h, 0) for h in ML); out['share_samples'] = [r['counts'].get(h, 0) / cm for h in ML]
    if noise == 'noisy':
        out['share_mitigated'] = shares(mitigate(r['exact'], name), ML)[0].tolist()
        tot = sum(r['counts'].values())
        out['share_samples_mitigated'] = shares(mitigate({h: v / tot for h, v in r['counts'].items()}, name), ML)[0].tolist()
    return out


for c in CASES:
    name = c['instance']; P = layout(name); w = c['weights']; ML = c['labels']; W = np.array(c['history_weight'], float)
    top = c['top']; R = res.setdefault(c['case'], {})
    for T in TS:
        for noise, read in READS:
            K = constants(T, noise, read)
            enc = {'raw': (raw_encoding(name, w), None, {})}
            bias = [K['zero'][name][h] for h in ML]
            for tag, theta, beta, use0 in (('th', c['theta_pred'], K['beta_th'], False), ('th0', c['theta_pred'], K['beta_th'], True),
                                           ('lai0', c['theta_pred'], K['beta_lai'], True), ('kl0', c['theta_kl'], K['beta_th'], True)):
                what, s_used, info = encode_qp(name, w, theta / beta, P, bias=bias if use0 else None, beta=beta)
                enc[tag] = (what, dict(theta=theta, beta=beta, s_req=theta / beta, s_used=s_used, theta_eff=s_used * beta,
                                       zero_point=use0), info)
            for tag, (what, meta, info) in enc.items():
                key = f'T{T}|{noise}|{read}|{tag}'
                if key in R: continue
                if tag == 'raw' and read == 'mit':     # same physical run as the raw-read one; only the read differs
                    R[key] = dict(R[f'T{T}|noisy|raw|raw']); R[key]['share_used'] = R[key]['share_mitigated']; continue
                t0 = time.time(); out = run(f"{c['case']}|{key}", name, P, what, T, noise)
                eps = energies(name, P, DELTAF * np.asarray(what)); gs = sorted(eps.items(), key=lambda t: -t[1])
                out.update(T=T, noise=noise, read=read, encoding=tag, what=list(map(float, what)), meta=meta,
                           delta=info.get('delta'), saturated=info.get('saturated'),
                           ground=dict(state=gs[0][0], is_top=gs[0][0] == top, gap=gs[0][1] - gs[1][1], next=gs[1][0]),
                           eps_maximal=[eps[h] for h in ML])
                out['share_used'] = out[SK[read]]
                R[key] = out
                print(c['case'], key, f'{time.time() - t0:.1f}s', 'boot', [round(x, 3) for x in c['bootstrap']], 'dev',
                      [round(x, 3) for x in out['share_used']], 'TV', round(tv(out['share_used'], c['bootstrap']), 3),
                      'sat' if out['saturated'] else '', flush=True)
    json.dump(res, open(outp, 'w'), indent=0)
print('done')
