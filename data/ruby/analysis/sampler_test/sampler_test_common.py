"""Shared definitions for the sampler_test_* scripts (23 Nov 2026): does a WEIGHTED neutral-atom device sample
stemmatic uncertainty? Emulation only, Mac Studio, no hardware. Read-only on every existing file.

Instances = the three Ruby registers (single copy c0), graphs copied from check_common.py (frozen scorers).
Drive = ruby_mis_prep.py / emu_dist.py exactly (Omega 2pi*1.5, delta -2pi*4 -> 2pi*0.9 interpolated, 4 knots;
DMM detuning map 1 - w_eff, ramp 0 -> -delta_f; DMM bottom -60, total -600), with T and the weight floor a varied.
Noise = the preps' vendor-nominal model: p_fp 0.02, p_fn 0.10, dephasing 1/4.5 us^-1.
Run with pulser 1.9.1 and qutip 5.3.1."""
import json, os, itertools as it, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'sampler_test_20261123')
os.makedirs(OUT, exist_ok=True)
C6 = 865723.02
TWO_PI = 2 * math.pi
OMEGA = TWO_PI * 1.5
DELTA0, DELTAF = -TWO_PI * 4.0, TWO_PI * 0.9
RB = (C6 / OMEGA) ** (1 / 6)

INST = {
    'priscus': dict(
        fam=['BMP', 'MP', 'EB', 'EP', 'EMP'],
        W={'BMP': 45, 'MP': 8, 'EB': 4, 'EP': 3, 'EMP': 3},
        conflicts=[('BMP', 'EB'), ('BMP', 'EP'), ('BMP', 'EMP'), ('MP', 'EP'), ('EB', 'EP'), ('EB', 'EMP')],
        prep='ruby_mis_prep.json'),
    'lai': dict(
        fam=['CG', 'DE', 'AB', 'DEF', 'EF', 'CEG', 'DF'],
        W={'CG': 176, 'DE': 43, 'AB': 38, 'DEF': 24, 'EF': 23, 'CEG': 13, 'DF': 13},
        conflicts=[('DE', 'EF'), ('DE', 'CEG'), ('DE', 'DF'), ('DEF', 'CEG'), ('EF', 'CEG'), ('EF', 'DF')],
        prep='lai_mis_prep.json'),
    'lawrence': dict(
        fam=['B+F+Mi+P+R', 'B+R', 'P+R', 'B+F+P+R', 'B+F+R'],
        W={'B+F+Mi+P+R': 11, 'B+R': 4, 'P+R': 4, 'B+F+P+R': 3, 'B+F+R': 3},
        conflicts=[('P+R', 'B+R'), ('P+R', 'B+F+R')],
        prep='lawrence_mis_prep.json'),
}


def conf(g):
    C = {frozenset(e) for e in g['conflicts']}
    return lambda a, b: frozenset((a, b)) in C


def independent_sets(g):
    c = conf(g); F = g['fam']
    return [frozenset(s) for k in range(len(F) + 1) for s in it.combinations(F, k)
            if all(not c(a, b) for a, b in it.combinations(s, 2))]


def maximal_sets(g):
    ind = independent_sets(g)
    return [s for s in ind if not any(s < t for t in ind)]


def lab(g, s):
    return '{' + ','.join(x for x in g['fam'] if x in s) + '}'


def weff(g, a, W=None):
    W = W or g['W']; m = max(W.values())
    return {q: a + (1 - a) * W[q] / m for q in g['fam']}


def register_positions(name):
    """c0 copy of the frozen register, centred (as in ruby_mis_prep.py: one copy = ext coordinates)."""
    g = INST[name]; reg = json.load(open(os.path.join(HERE, g['prep'])))['register_um']
    P = np.array([reg[f'c0_{f}'] for f in g['fam']], float)
    return P - P.mean(0)


def classical_spectrum(g, P, delta_vec):
    """End-of-sweep classical Rydberg energy (Omega = 0) with ALL vdW tails, over all 2^n configurations."""
    F = g['fam']; n = len(F)
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf); V = C6 / D ** 6
    Es = []
    for bits in it.product([0, 1], repeat=n):
        z = np.array(bits); Es.append((float(0.5 * z @ V @ z - delta_vec @ z), frozenset(F[i] for i in range(n) if bits[i])))
    Es.sort(key=lambda t: t[0])
    return Es


def ground_check(g, P, a, W=None):
    we = weff(g, a, W); dv = np.array([DELTAF * we[f] for f in g['fam']])
    Es = classical_spectrum(g, P, dv)
    ind = independent_sets(g); mwis = max(ind, key=lambda s: sum((W or g['W'])[x] for x in s))
    return dict(ground=lab(g, Es[0][1]), is_mwis=Es[0][1] == mwis, gap=Es[1][0] - Es[0][0],
                next=lab(g, Es[1][1]), mwis=lab(g, mwis))


def geometry(g, P):
    c = conf(g); F = g['fam']; n = len(F)
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1))
    dc = [D[i, j] for i, j in it.combinations(range(n), 2) if c(F[i], F[j])]
    dn = [D[i, j] for i, j in it.combinations(range(n), 2) if not c(F[i], F[j])]
    return dict(min_spacing=float(min(D[i, j] for i, j in it.combinations(range(n), 2))),
                max_conflict=float(max(dc)), min_compatible=float(min(dn)),
                conflict_over_Rb=float(max(dc) / RB), compatible_over_Rb=float(min(dn) / RB),
                V_conflict_min=float(C6 / max(dc) ** 6), V_compatible_max=float(C6 / min(dn) ** 6),
                radius=float(np.sqrt((P ** 2).sum(1)).max()))
