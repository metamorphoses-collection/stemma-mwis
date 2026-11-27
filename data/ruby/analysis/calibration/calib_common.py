"""Shared helpers for the calib_* scripts (23 Nov 2026): can a weighted neutral-atom device be CALIBRATED so that its
outcome distribution over stemmatic histories reproduces the bootstrap uncertainty of the evidence?
Emulation only (Mac Studio). Reuses sampler_test_common.py (instances, conflict graphs, layouts, classical spectrum).
Read-only on every existing file; outputs in calib_20261123/.
Run with pulser 1.9.1 and qutip 5.3.1."""
import json, os, math, itertools as it, dataclasses, warnings
import numpy as np
from scipy.stats import norm, skellam
from sampler_test_common import (INST, HERE, OUT as SOUT, C6, OMEGA, DELTA0, DELTAF, RB, conf, independent_sets, maximal_sets,
                                 lab, register_positions, classical_spectrum, geometry)

warnings.filterwarnings('ignore')
OUT = os.path.join(HERE, 'calib_20261123'); os.makedirs(OUT, exist_ok=True)
STEM = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'lawrence'))
LAI_UNITS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'lai_ombre', 'units_v2.json')
LOGIT_PROBIT = 1.702   # logistic(1.702 z) ~ Phi(z), max abs error 0.0095
FP, FN, DEPH = 0.02, 0.10, 1 / 4.5


# ------------------------------------------------------------------ units (support vectors per unit, register families only)
def units_priscus():
    tally = {'BMP': 45, 'MP': 8, 'EB': 4, 'EP': 3, 'EMP': 3}
    return [{f} for f, k in tally.items() for _ in range(k)] + [set(), set()]


def units_lai():
    """units_v2.json, unrooted split supports restricted to the 7 register families (sampler_test_bootstrap.py logic)."""
    ALL = frozenset('ABCDEFG')
    key = {'CG': 'CG', 'DE': 'DE', 'AB': 'CDEFG', 'DEF': 'DEF', 'EF': 'EF', 'CEG': 'CEG', 'DF': 'DF'}
    out = []
    for u in json.load(open(LAI_UNITS)):
        s = set()
        for gr in u['groups'][1:]:
            if gr and len(gr) >= 2 and frozenset(gr) != ALL:
                g = frozenset(gr)
                if 2 <= len(g) <= 5: s.add(g if 'A' not in g else ALL - g)
        out.append({r for r, k in key.items() if frozenset(k) in s})
    return out


def units_lawrence():
    L = [l for l in open(os.path.join(STEM, 'lawrence_bhl4759.num')).read().split('\n') if l.strip()]
    wit = L[0].split(); rows = [l.split() for l in L[1:]]
    G = json.load(open(os.path.join(STEM, 'lawrence_native_graph_s3.json')))
    famsets = [frozenset(x.split('+')) for x in G['names']]
    out = []
    for r in rows:
        obs = [i for i in range(len(wit)) if r[i] != '?']; s = set()
        for v in set(r[i] for i in obs):
            A = frozenset(wit[i] for i in obs if r[i] == v)
            if 'T' in A: A = frozenset(wit[i] for i in obs) - A
            if A in famsets: s.add(G['names'][famsets.index(A)])
        out.append(s & set(INST['lawrence']['fam']))
    return out


UNITS = {'priscus': units_priscus, 'lai': units_lai, 'lawrence': units_lawrence}


def X_of(name, U):
    F = INST[name]['fam']; return np.array([[f in u for f in F] for u in U], dtype=np.int64)


def M_of(name):
    g = INST[name]; mx = maximal_sets(g)
    return [lab(g, s) for s in mx], np.array([[f in s for f in g['fam']] for s in mx], dtype=np.int64)


# ------------------------------------------------------------------ bootstrap and its analytic approximations
def bootstrap(X, M, B=20000, seed=1):
    """Share of each maximal history being the max-weight one over B unit resamplings (ties split equally)."""
    rng = np.random.default_rng(seed); n = len(X); P = np.zeros(len(M))
    C = rng.multinomial(n, np.full(n, 1 / n), size=B)            # B x n resample counts
    V = (C @ X) @ M.T                                             # B x H history scores
    mx = V.max(1, keepdims=True); top = (V == mx); P = (top / top.sum(1, keepdims=True)).sum(0) / B
    return P


def pair_stats(X, M, h, k):
    """D = v_h - v_k on the data; exact multinomial bootstrap variance sigma^2 = sum d_u^2 - D^2/N; N_diff; max|d_u|."""
    d = X @ (M[h] - M[k]); N = len(d); D = int(d.sum())
    return dict(D=D, var=float((d ** 2).sum() - D ** 2 / N), N_diff=int((d != 0).sum()), max_abs_d=int(np.abs(d).max()),
                n_plus=int((d > 0).sum()), n_minus=int((d < 0).sum()))


def p_pair_normal(st):
    return float(norm.cdf(st['D'] / math.sqrt(st['var']))) if st['var'] > 0 else (1.0 if st['D'] > 0 else 0.5 if st['D'] == 0 else 0.0)


def p_pair_skellam(st):
    """Poisson limit (exclusive units only): K+ ~ Pois(n+), K- ~ Pois(n-); P(K+ > K-) + P(=)/2."""
    a, b = st['n_plus'], st['n_minus']
    if a == 0 and b == 0: return 0.5
    if b == 0: return 1 - 0.5 * math.exp(-a)
    if a == 0: return 0.5 * math.exp(-b)
    return float(skellam.sf(0, a, b) + 0.5 * skellam.pmf(0, a, b))


def mvn_probit(X, M, B=200000, seed=2):
    """Gaussian (multinomial-probit) approximation of the bootstrap: v* ~ N(M w, M Sigma M^T), Sigma = X^T X - w w^T / N."""
    w = X.sum(0); N = len(X); S = X.T @ X - np.outer(w, w) / N
    mu = M @ w; C = M @ S @ M.T
    rng = np.random.default_rng(seed); V = rng.multivariate_normal(mu, C, size=B, method='eigh')
    return np.bincount(V.argmax(1), minlength=len(M)) / B


def softmax(s):
    s = np.asarray(s, float); e = np.exp(s - s.max()); return e / e.sum()


def tv(p, q):
    return 0.5 * float(np.abs(np.asarray(p, float) - np.asarray(q, float)).sum())


# ------------------------------------------------------------------ layouts
def layout(name):
    """The checked layouts of the sampler test: Priscus frozen register c0; Lai opt20.0; Lawrence opt12.0."""
    F = INST[name]['fam']
    if name == 'priscus': return register_positions('priscus')
    key = {'lai': '20.0', 'lawrence': '12.0'}[name]
    L = json.load(open(os.path.join(SOUT, f'layout_{name}.json')))['optimised'][key]
    return np.array([L['positions_um'][f] for f in F], float)


def interactions(P):
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf); return C6 / D ** 6


def energies(name, P, delta):
    """Classical end-of-sweep 'reward' eps_S = sum_{i in S} delta_i - sum_{i<j in S} V_ij for every bitstring label."""
    g = INST[name]; F = g['fam']; n = len(F); V = interactions(P); out = {}
    for bits in it.product([0, 1], repeat=n):
        z = np.array(bits, float); out[lab(g, {F[i] for i in range(n) if bits[i]})] = float(delta @ z - 0.5 * z @ V @ z)
    return out


def ground_check_delta(name, P, delta, target_label):
    """Lowest classical energy state (max reward) with all tails, vs the target (max-weight) history."""
    e = energies(name, P, delta); s = sorted(e.items(), key=lambda t: -t[1])
    return dict(ground=s[0][0], is_target=s[0][0] == target_label, gap=s[0][1] - s[1][1], next=s[1][0])


# ------------------------------------------------------------------ emulation with an arbitrary per-atom final detuning
def make_device(off_device=False):
    from pulser.devices import AnalogDevice
    from pulser.channels.dmm import DMM
    Dev = dataclasses.replace(AnalogDevice.to_virtual(), name='RubyLike+DMM', max_radial_distance=50, max_atom_num=100,
                              min_atom_distance=5, dmm_objects=(DMM(bottom_detuning=-60, total_bottom_detuning=-600),))
    if off_device: Dev = dataclasses.replace(Dev, max_sequence_duration=None)
    return Dev


def emulate(name, P, what, T, noisy, deltaf=DELTAF, samples=4000, seed=0):
    """what: None (unweighted, every atom ends at deltaf) or a vector w_hat in [0,1] (final local detuning deltaf * w_hat).
    Same drive as sampler_test_emulate.py (Omega 4 knots; delta d0, d0/2, 0.6 df, df; DMM map 1 - w_hat ramped 0 -> -df).
    Returns exact outcome distribution (post readout channel if noisy), 4000-sample counts, both keyed by history label."""
    from pulser import Register, Sequence, Pulse, NoiseModel
    from pulser.waveforms import InterpolatedWaveform, RampWaveform
    from pulser_simulation import QutipEmulator, SimConfig
    g = INST[name]; F = g['fam']; n = len(F)
    reg = Register({f: tuple(map(float, P[i])) for i, f in enumerate(F)})
    s = Sequence(reg, make_device(T > 6000)); s.declare_channel('ryd', 'rydberg_global')
    if what is not None:
        what = np.clip(np.asarray(what, float), 0, 1)
        s.config_detuning_map(reg.define_detuning_map({f: float(1 - what[i]) for i, f in enumerate(F)}, 'dmap'), 'dmm_0')
    s.add(Pulse(InterpolatedWaveform(T, [1e-9, OMEGA, OMEGA, 1e-9]),
                InterpolatedWaveform(T, [DELTA0, DELTA0 * 0.5, deltaf * 0.6, deltaf]), 0), 'ryd')
    if what is not None and np.any(what < 1):
        s.add_dmm_detuning(RampWaveform(T, 0, -deltaf), 'dmm_0')
    vendor = NoiseModel(p_false_pos=FP, p_false_neg=FN, dephasing_rate=DEPH, runs=40, samples_per_run=50)
    sim = QutipEmulator.from_sequence(s, sampling_rate=0.05, config=SimConfig.from_noise_model(vendor) if noisy else None)
    r = sim.run()
    np.random.seed(seed); cnt = r.sample_final_state(N_samples=samples)
    M = r.get_final_state().full()
    p = np.real(np.diag(M)) if M.shape[0] == M.shape[1] else np.abs(M[:, 0]) ** 2
    BITS = [''.join(b) for b in it.product('01', repeat=n)]
    pre = {}
    for k, v in enumerate(p):
        pre[''.join('1' if ch == '0' else '0' for ch in format(k, f'0{n}b'))] = float(v)
    ex = readout(pre, n) if noisy else pre
    L = lambda b: lab(g, {F[i] for i in range(n) if b[i] == '1'})
    return dict(exact={L(b): v for b, v in ex.items()}, pre_readout={L(b): v for b, v in pre.items()},
                counts={L(b): int(v) for b, v in cnt.items()})


def confusion(n):
    m = np.array([[1 - FP, FN], [FP, 1 - FN]])   # rows: observed 0/1, cols: true 0/1
    R = np.array([[1.0]])
    for _ in range(n): R = np.kron(R, m)
    return R


def readout(pre, n):
    BITS = [''.join(b) for b in it.product('01', repeat=n)]
    v = np.array([pre.get(b, 0.0) for b in BITS]); o = confusion(n) @ v
    return {b: float(x) for b, x in zip(BITS, o)}


def mitigate(dist_by_label, name):
    """Invert the (known, vendor-nominal) independent readout channel; clip negatives, renormalise.
    Works on an exact distribution or on empirical frequencies."""
    g = INST[name]; F = g['fam']; n = len(F)
    BITS = [''.join(b) for b in it.product('01', repeat=n)]
    L = lambda b: lab(g, {F[i] for i in range(n) if b[i] == '1'})
    v = np.array([dist_by_label.get(L(b), 0.0) for b in BITS], float); v = v / v.sum()
    t = np.linalg.solve(confusion(n), v); t = np.clip(t, 0, None); t /= t.sum()
    return {L(b): float(x) for b, x in zip(BITS, t)}


def shares(dist, ML):
    p = np.array([dist.get(h, 0.0) for h in ML], float); s = p.sum()
    return (p / s if s > 0 else np.full(len(ML), np.nan)), float(s)
