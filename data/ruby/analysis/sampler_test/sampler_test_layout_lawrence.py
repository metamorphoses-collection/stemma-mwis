"""Task 2: a 5-atom Lawrence-core layout whose weighted end-of-sweep lowest-energy state (Omega = 0, ALL vdW tails) is the
max-weight set E = {B+F+Mi+P+R, B+R, B+F+P+R, B+F+R}, with a margin, under the preps' constraints:
min spacing >= 5 um; conflicting pairs <= 0.90 R_b; compatible pairs >= 1.25 R_b; R_b = (C6/Omega)^(1/6) = 6.72 um.
Objective: maximise the smaller of the two weighted gaps (floor a = 0.3 and a = 0.1), for a set of footprint radii
(a smaller footprint tiles more copies). Random multistart + Nelder-Mead on a penalised objective; constraints and
the ground state are then re-checked exactly over all 32 configurations.
Also tabulates the same ground-state check for the frozen Priscus and Lai registers at a = 0.3 and 0.1.
Output: sampler_test_20261123/layout_lawrence.json"""
import json, os, itertools as it
import numpy as np
from scipy.optimize import minimize
from sampler_test_common import INST, OUT, RB, C6, DELTAF, conf, register_positions, ground_check, geometry, weff, classical_spectrum, lab

g = INST['lawrence']; F = g['fam']; n = 5; c = conf(g)
pairs = list(it.combinations(range(n), 2)); isc = [c(F[i], F[j]) for i, j in pairs]
A_LIST = (0.3, 0.1)


Z = np.array(list(it.product([0, 1], repeat=n)), float)
E_IDX = int(np.flatnonzero((Z == np.array([1, 1, 0, 1, 1])).all(1))[0])   # the max-weight set E
DV = {a: np.array([DELTAF * weff(g, a)[f] for f in F]) for a in A_LIST}


def gaps(P):
    """Fast vectorised version of ground_check (re-checked exactly with ground_check afterwards):
    min over floors of E(next best) - E(E); negative if E is not the lowest."""
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf); V = C6 / D ** 6
    tail = 0.5 * np.einsum('ki,ij,kj->k', Z, V, Z)
    out = []
    for a in A_LIST:
        Eng = tail - Z @ DV[a]; e0 = Eng[E_IDX]; others = np.delete(Eng, E_IDX)
        out.append(others.min() - e0)
    return min(out)


def penalty(P, rfp):
    pen = 0.0
    for (i, j), cf in zip(pairs, isc):
        d = np.linalg.norm(P[i] - P[j]); pen += max(0, 5.0 - d) ** 2
        pen += max(0, d - 0.90 * RB) ** 2 if cf else max(0, 1.25 * RB - d) ** 2
    Q = P - P.mean(0); pen += max(0, np.sqrt((Q ** 2).sum(1)).max() - rfp) ** 2
    return pen


def obj(x, rfp):
    P = x.reshape(n, 2); return -gaps(P) + 50 * penalty(P, rfp)


rng = np.random.default_rng(20261123)
out = {'constraints': dict(min_spacing_um=5.0, conflict_max_over_Rb=0.90, compatible_min_over_Rb=1.25, R_b_um=RB,
                           delta_f=DELTAF, weight_floors=A_LIST)}
P0 = register_positions('lawrence')
out['frozen_register'] = dict(geometry=geometry(g, P0), **{f'a{a}': ground_check(g, P0, a) for a in A_LIST})
best_by_r = {}
for rfp in (8.0, 10.0, 12.0, 15.0, 20.0):
    best = None
    for k in range(80):
        x0 = rng.uniform(-rfp, rfp, 2 * n)
        r = minimize(obj, x0, args=(rfp,), method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-4, fatol=1e-6))
        P = r.x.reshape(n, 2); P = P - P.mean(0)
        if penalty(P, rfp + 1e-3) > 1e-6: continue   # feasible to ~1e-3 um
        gp = gaps(P)
        if best is None or gp > best[0]: best = (gp, P)
    if best is None:
        best_by_r[str(rfp)] = dict(footprint_radius_um=rfp, feasible=False); print(rfp, 'no feasible layout found', flush=True); continue
    gp, P = best
    rec = dict(footprint_radius_um=rfp, min_gap=gp, positions_um={F[i]: [round(float(v), 3) for v in P[i]] for i in range(n)},
               geometry=geometry(g, P), unweighted=None)
    for a in A_LIST: rec[f'a{a}'] = ground_check(g, P, a)
    Es = classical_spectrum(g, P, np.full(n, DELTAF)); rec['unweighted'] = dict(ground=lab(g, Es[0][1]), gap=Es[1][0] - Es[0][0])
    best_by_r[str(rfp)] = rec
    print(rfp, round(gp, 3), rec['geometry'], {a: (rec[f'a{a}']['ground'], round(rec[f'a{a}']['gap'], 3)) for a in A_LIST}, flush=True)
out['optimised'] = best_by_r
# upper bound: the gap can never exceed the smallest local detuning of E's atoms (removing it costs delta_f * w_eff_min)
out['gap_upper_bound'] = {f'a{a}': DELTAF * min(weff(g, a).values()) for a in A_LIST}

# other registers, both floors
for name in ('priscus', 'lai'):
    P = register_positions(name)
    out[f'{name}_frozen_register'] = dict(geometry=geometry(INST[name], P), **{f'a{a}': ground_check(INST[name], P, a) for a in A_LIST})
json.dump(out, open(os.path.join(OUT, 'layout_lawrence.json'), 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in out.items() if k != 'optimised'}, indent=1, default=float))
