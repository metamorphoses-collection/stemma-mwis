"""Extra (not in the brief): the frozen Lai register's weighted ground state is the MWIS only at a = 0.3, by 0.03 rad/us,
and drops AB at a = 0.1 / 0.0 (tails). Simply moving the isolated CG and AB away is not enough (DEF then drops at
a = 0.1: its local detuning 0.223*delta_f = 1.26 rad/us is below its tails). So: optimise all 7 positions as
sampler_test_layout_lawrence.py does (min spacing 5 um, conflicts <= 0.90 R_b, compatible >= 1.25 R_b, footprint
radius <= R), maximising the smaller weighted gap over a in {0.3, 0.1}; ground state re-checked exactly (128 configs).
Output: sampler_test_20261123/layout_lai.json"""
import json, os, itertools as it
import numpy as np
from scipy.optimize import minimize
from sampler_test_common import INST, OUT, RB, C6, DELTAF, conf, register_positions, ground_check, geometry, weff, independent_sets

g = INST['lai']; F = g['fam']; n = len(F); c = conf(g)
pairs = list(it.combinations(range(n), 2)); isc = [c(F[i], F[j]) for i, j in pairs]
A_LIST = (0.3, 0.1)
Z = np.array(list(it.product([0, 1], repeat=n)), float)
mw = max(independent_sets(g), key=lambda s: sum(g['W'][x] for x in s))
E_IDX = int(np.flatnonzero((Z == np.array([f in mw for f in F], float)).all(1))[0])
DV = {a: np.array([DELTAF * weff(g, a)[f] for f in F]) for a in A_LIST}


def gaps(P):
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf); V = C6 / D ** 6
    tail = 0.5 * np.einsum('ki,ij,kj->k', Z, V, Z)
    return min((np.delete(tail - Z @ DV[a], E_IDX).min() - (tail - Z @ DV[a])[E_IDX]) for a in A_LIST)


def penalty(P, rfp):
    pen = 0.0
    for (i, j), cf in zip(pairs, isc):
        d = np.linalg.norm(P[i] - P[j]); pen += max(0, 5.0 - d) ** 2
        pen += max(0, d - 0.90 * RB) ** 2 if cf else max(0, 1.25 * RB - d) ** 2
    Q = P - P.mean(0); pen += max(0, np.sqrt((Q ** 2).sum(1)).max() - rfp) ** 2
    return pen


def obj(x, rfp):
    P = x.reshape(n, 2); return -gaps(P) + 1e4 * penalty(P, rfp)


rng = np.random.default_rng(20261123)
# hand-built feasible seed (random starts in 14 dims never reached feasibility): DE/EF a conflicting pair, CEG and DF on
# either side of it (both conflict with DE and EF, compatible with each other), DEF beyond CEG, CG and AB far away.
SEED = np.array([[0, 15], [0, 2.6], [0, -15], [-10.8, 0], [0, -2.6], [-5.3, 0], [5.3, 0]], float)  # order = F
assert [*F] == ['CG', 'DE', 'AB', 'DEF', 'EF', 'CEG', 'DF']
out = {'gap_upper_bound': {f'a{a}': DELTAF * min(weff(g, a).values()) for a in A_LIST},
       'frozen_register': {f'a{a}': ground_check(g, register_positions('lai'), a) for a in (0.3, 0.1, 0.0)}, 'optimised': {}}
for rfp in (16.0, 20.0):
    best = None
    for k in range(80):
        r = minimize(obj, (SEED + rng.normal(0, 0.7 if k else 0.0, SEED.shape)).ravel(), args=(rfp,), method='Nelder-Mead',
                     options=dict(maxiter=20000, xatol=1e-4, fatol=1e-7))
        P = r.x.reshape(n, 2); P = P - P.mean(0)
        if penalty(P, rfp + 1e-3) > 1e-6: continue
        gp = gaps(P)
        if best is None or gp > best[0]: best = (gp, P)
    if best is None:
        out['optimised'][str(rfp)] = dict(feasible=False); print(rfp, 'none', flush=True); continue
    gp, P = best
    ms = min(np.linalg.norm(P[i] - P[j]) for i, j in pairs)
    if ms < 5.0001: P = P * (5.0001 / ms); gp = gaps(P)   # NM sits on the 5 um bound to ~1e-5: rescale to clear Pulser's 1e-6 check
    rec = dict(footprint_radius_um=rfp, min_gap=gp, positions_um={f: P[i].tolist() for i, f in enumerate(F)}, geometry=geometry(g, P),
               ground_checks={f'a{a}': ground_check(g, P, a) for a in (0.3, 0.1, 0.0)})
    out['optimised'][str(rfp)] = rec
    print(rfp, round(gp, 3), rec['geometry'], {a: (v['ground'], round(v['gap'], 3)) for a, v in rec['ground_checks'].items()}, flush=True)
json.dump(out, open(os.path.join(OUT, 'layout_lai.json'), 'w'), indent=1, default=float)
print(json.dumps(out['gap_upper_bound']))
