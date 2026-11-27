"""Test 2: graded instances with known bootstrap uncertainty.
Cases = unit subsamples of Lawrence (75 units) and Lai (units_v2.json, 843 units, restricted to the 7 register families),
chosen so that the bootstrap share of the top history spans ~0.55 .. 0.99, plus Priscus as is (and the two full sets).
Conflict graphs (and so the maximal histories) and layouts are those of the registers (layout_lawrence opt12.0,
layout_lai opt20.0, Priscus register c0); only the family weights change with the subsample.
Selection: random subsets (seeded) of several sizes; for each target share the subset whose quick bootstrap (B=4000)
is closest is kept, then re-bootstrapped with B=20000. Cases with a tie for the top weight are excluded.
Per case: weights, bootstrap shares, exact pair statistics (D, sigma_D^2 = sum d_u^2 - D^2/L, N_diff, max|d_u|) for the top
history vs every other, normal (Phi(D/sigma)) and multinomial-probit (MVN) approximations, and the softmax-encoding
ceilings: theta_pred = 1.702/sigma_D(top vs runner-up) (no bootstrap used), theta_KL = argmin KL(p_boot||softmax(theta W)),
theta_TV = argmin TV. Output: calib_20261123/cases.json"""
import json, os
import numpy as np
from scipy.optimize import minimize_scalar
from calib_common import (INST, OUT, UNITS, X_of, M_of, bootstrap, pair_stats, p_pair_normal, mvn_probit, softmax, tv, LOGIT_PROBIT)

rng = np.random.default_rng(20261123)
TARGETS = [0.55, 0.62, 0.70, 0.78, 0.86, 0.93, 0.99]


def describe(name, idx, U, tag):
    X = X_of(name, [U[i] for i in idx]); ML, M = M_of(name); w = X.sum(0); v = M @ w
    o = np.argsort(-v, kind='stable'); h = int(o[0])
    P = bootstrap(X, M, 20000, 1)
    pairs = {}
    for k in o[1:]:
        st = pair_stats(X, M, h, int(k)); st['p_normal'] = p_pair_normal(st)
        # two-history bootstrap restricted to the pair (for the pairwise check)
        Pk = bootstrap(X, M[[h, int(k)]], 20000, 3); st['p_boot_pair'] = float(Pk[0]); pairs[ML[int(k)]] = st
    runner = ML[int(o[1])]; sig = np.sqrt(pairs[runner]['var'])
    th_pred = LOGIT_PROBIT / sig
    def kl(t):
        q = softmax(t * v); return float(sum(p * np.log(p / q[i]) for i, p in enumerate(P) if p > 0))
    th_kl = float(minimize_scalar(kl, bounds=(1e-4, 20), method='bounded').x)
    th_tv = float(minimize_scalar(lambda t: tv(P, softmax(t * v)), bounds=(1e-4, 20), method='bounded').x)
    return dict(case=tag, instance=name, units=len(idx), unit_index=[int(i) for i in idx],
                weights=dict(zip(INST[name]['fam'], w.tolist())), labels=ML, history_weight=v.tolist(), top=ML[h],
                runner_up=runner, bootstrap=P.tolist(), mvn=mvn_probit(X, M).tolist(), pairs=pairs,
                sigma_top_runner=float(sig), z_top_runner=float(pairs[runner]['D'] / sig),
                theta_pred=th_pred, theta_kl=th_kl, theta_tv=th_tv,
                softmax_pred=softmax(th_pred * v).tolist(), softmax_kl=softmax(th_kl * v).tolist(),
                tv_softmax_pred=tv(P, softmax(th_pred * v)), tv_softmax_kl=tv(P, softmax(th_kl * v)), tv_mvn=tv(P, mvn_probit(X, M)))


cases = []
for name, sizes in (('lawrence', [20, 25, 30, 35, 40, 50, 60, 70]), ('lai', [30, 50, 80, 120, 180, 260, 400, 600])):
    U = UNITS[name](); X = X_of(name, U); ML, M = M_of(name); pool = []
    for n in sizes:
        for _ in range(60):
            idx = np.sort(rng.choice(len(U), n, replace=False)); v = M @ X[idx].sum(0); s = np.sort(v)[::-1]
            if s[0] == s[1]: continue
            pool.append((idx, bootstrap(X[idx], M, 4000, 5).max()))
    used = set()
    for t in TARGETS:
        best = min((p for p in pool if id(p[0]) not in used), key=lambda p: abs(p[1] - t))
        if abs(best[1] - t) > 0.04: print(name, 'target', t, 'nearest', round(best[1], 3), '(skipped)'); continue
        used.add(id(best[0])); cases.append(describe(name, best[0], U, f'{name}_sub{len(best[0])}_t{t:.2f}'))
    cases.append(describe(name, np.arange(len(U)), U, f'{name}_full'))
U = UNITS['priscus'](); cases.append(describe('priscus', np.arange(len(U)), U, 'priscus_full'))

json.dump(cases, open(os.path.join(OUT, 'cases.json'), 'w'), indent=1)
for c in cases:
    print(f"{c['case']:28s} n={c['units']:4d} top={c['top']:34s} boot_top={max(c['bootstrap']):.3f} z={c['z_top_runner']:.2f} "
          f"theta_pred={c['theta_pred']:.3f} theta_KL={c['theta_kl']:.3f} TV(softmax_pred)={c['tv_softmax_pred']:.3f} "
          f"TV(softmax_KL)={c['tv_softmax_kl']:.3f} TV(mvn)={c['tv_mvn']:.3f}")
    for k, st in c['pairs'].items():
        print(f"     vs {k:34s} D={st['D']:3d} var={st['var']:.2f} Ndiff={st['N_diff']} max|d|={st['max_abs_d']} "
              f"normal={st['p_normal']:.3f} boot_pair={st['p_boot_pair']:.3f}")
