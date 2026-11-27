"""DISTRIBUTION variant of ruby_mis_prep.py (same register, drive and noise model): keeps the full measured counts and
scans the sweep duration T (env SWEEP_NS, comma list). Emulation only. Original docstring follows.
Priscus stemma MIS on Ruby (no addressing) + MWIS emulation at Ruby geometry. PREPARATION ONLY: builds,
checks and emulates; submits nothing.

Device limits asserted below (neutral-atom QPU "Ruby", Oct-Nov 2026):
C6 = 865723.02 rad/us um^6 (ANGULAR, n=60) · Omega <= 15.708 rad/us (QPU-enforced) · |delta| <= 62.83 rad/us
· duration <= 6000 ns · min spacing 5 um · radius <= 50 um · <= 100 atoms
· min_avg_amp 0.3142 · margin: compatible pairs >= 1.25 R_b (R_b at the drive Omega).

usage: python3 ruby_mis_prep.py <layout6.json> <outdir>
"""
import json, sys, itertools as it, dataclasses
import numpy as np
from pulser import Register, Sequence, Pulse
from pulser.waveforms import InterpolatedWaveform, RampWaveform
from pulser.devices import AnalogDevice
from pulser.channels.dmm import DMM
from pulser_simulation import QutipEmulator, SimConfig
from pulser import NoiseModel

lay, OUT = sys.argv[1], sys.argv[2]
L = json.load(open(lay)); names, W, X, adj = L['names'], L['W'], np.array(L['X']), np.array(L['adj'])
n = len(names); iu = np.triu_indices(n, 1)
C6 = 865723.02
assert abs(AnalogDevice.interaction_coeff - C6) < 1, AnalogDevice.interaction_coeff

# ---- drive and geometry (frozen) ----
TWO_PI = 2 * np.pi
OMEGA = TWO_PI * 1.5            # 9.42 rad/us  (cap 15.708)
DELTA0, DELTAF = -TWO_PI * 4.0, TWO_PI * 0.9
import os
TS = [int(x) for x in os.environ.get('SWEEP_NS', '6000').split(',')]
T = 6000                        # ns (preflight at the frozen 6000 ns; the scan below varies it)
SCALE = 7.2                     # um per layout unit (layout_ruby.json: min distance 0.72 -> 5.18 um)
A_COMPRESS = 0.3                # MWIS emulation: w_eff = a + (1-a) w/w_max (optimum asserted unchanged)
RB = (C6 / OMEGA) ** (1 / 6)
P = X * SCALE
D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1))
dconf, dcomp = D[iu][adj[iu] == 1], D[iu][adj[iu] == 0]
V = lambda d: C6 / d ** 6
PRE = {
    'Omega_rad_us': OMEGA, 'delta0': DELTA0, 'deltaf': DELTAF, 'T_ns': T, 'R_b_um': RB,
    'min_spacing_um': D[iu].min(), 'max_conflict_dist_um': dconf.max(), 'min_compatible_dist_um': dcomp.min(),
    'V_conflict_min': V(dconf.max()), 'V_compatible_max': V(dcomp.min()),
    'conflict_over_Rb': dconf.max() / RB, 'compatible_over_Rb': dcomp.min() / RB,
}
# ---- PREFLIGHT (must be able to fail) ----
assert OMEGA <= 15.708 and OMEGA >= 0.3142, 'Omega outside Ruby QPU limits'
assert max(abs(DELTA0), abs(DELTAF)) <= 62.83, 'detuning above Ruby cap'
assert 16 <= T <= 6000 and T % 4 == 0, 'duration'
assert PRE['min_spacing_um'] >= 5.0, 'atoms closer than 5 um'
assert PRE['conflict_over_Rb'] < 1.0, 'a conflicting pair is outside the blockade radius'
assert PRE['compatible_over_Rb'] >= 1.25, 'house margin r_nnn/R_b >= 1.25 violated'
assert PRE['V_conflict_min'] > 2 * DELTAF, 'final detuning not well inside the blockade'
assert PRE['V_compatible_max'] < DELTAF / 4, 'residual tail too strong for the final detuning'

# ---- MWIS weight compression for the DMM (emulation only) ----
WMAX = max(W.values()); WEFF = {q: A_COMPRESS + (1 - A_COMPRESS) * W[q] / WMAX for q in names}
_ind = [s_ for k_ in range(n + 1) for s_ in it.combinations(range(n), k_) if all(not adj[a, b] for a, b in it.combinations(s_, 2))]
_b0 = max(_ind, key=lambda s_: sum(W[names[i]] for i in s_)); _b1 = max(_ind, key=lambda s_: sum(WEFF[names[i]] for i in s_))
assert set(_b0) == set(_b1), 'weight compression changed the MWIS'
assert min(WEFF.values()) * DELTAF > PRE['V_compatible_max'], 'smallest local detuning below the residual tail'
PRE['WEFF'] = WEFF

# ---- tile copies on a ring inside the 50 um radius, >= 12 um atom-to-atom between copies ----
ext = P - P.mean(0); rad = np.sqrt((ext ** 2).sum(1)).max()
GAP = 12.0; ring = 50.0 - rad - 0.5
k = int(np.floor(np.pi / np.arcsin(min(1.0, (2 * rad + GAP) / (2 * ring)))))
centres = [(ring * np.cos(2 * np.pi * j / k), ring * np.sin(2 * np.pi * j / k)) for j in range(k)]
coords = {}
for k, (cx, cy) in enumerate(centres):
    for i in range(n):
        coords[f'c{k}_{names[i]}'] = (float(ext[i, 0] + cx), float(ext[i, 1] + cy))
allp = np.array(list(coords.values()))
dd = np.sqrt(((allp[:, None] - allp[None]) ** 2).sum(-1)); np.fill_diagonal(dd, 1e9)
assert dd.min() >= 5.0 and np.sqrt((allp ** 2).sum(1)).max() <= 50.0 and len(coords) <= 100
inter = min(np.linalg.norm(np.array(coords[f'c{a}_{names[i]}']) - np.array(coords[f'c{b}_{names[j]}']))
            for a, b in it.combinations(range(len(centres)), 2) for i in range(n) for j in range(n)) if len(centres) > 1 else None
PRE.update({'copies': len(centres), 'atoms': len(coords), 'min_intercopy_um': inter, 'V_intercopy_max': V(inter) if inter else None})

# ---- sequence (single copy for emulation; the tiled register for the run) ----
def seq_for(reg, device, weighted=False):
    s = Sequence(reg, device); s.declare_channel('ryd', 'rydberg_global')
    if weighted:
        dmap = reg.define_detuning_map({q: 1 - WEFF[q.split('_')[-1]] for q in reg.qubit_ids}, 'dmap')
        s.config_detuning_map(dmap, 'dmm_0')
    s.add(Pulse(InterpolatedWaveform(T, [1e-9, OMEGA, OMEGA, 1e-9]),
                InterpolatedWaveform(T, [DELTA0, DELTA0 * 0.5, DELTAF * 0.6, DELTAF]), 0), 'ryd')
    if weighted:
        s.add_dmm_detuning(RampWaveform(T, 0, -DELTAF), 'dmm_0')
    return s
RubyLike = dataclasses.replace(AnalogDevice.to_virtual(), name='RubyLike', max_radial_distance=50, max_atom_num=100, min_atom_distance=5)
tiled = Register(coords)
seq_tiled = seq_for(tiled, RubyLike)                # n=60 C6, 5 um, 50 um radius, global channel only
one = Register({f'c0_{names[i]}': tuple(ext[i]) for i in range(n)})

def stats(cnt):
    tot = sum(cnt.values())
    ok = lambda k: all(not (adj[a, b] and k[a] == '1' and k[b] == '1') for a, b in it.combinations(range(n), 2))
    ind = [s for k in range(n + 1) for s in it.combinations(range(n), k) if all(not adj[a, b] for a, b in it.combinations(s, 2))]
    m = max(len(s) for s in ind); MIS = [''.join('1' if i in s else '0' for i in range(n)) for s in ind if len(s) == m]
    best = max(ind, key=lambda s: sum(W[names[i]] for i in s)); bb = ''.join('1' if i in best else '0' for i in range(n))
    lab = lambda b: '{' + ','.join(names[i] for i in range(n) if b[i] == '1') + '}'
    return {'P_valid': sum(v for k, v in cnt.items() if ok(k)) / tot, 'P_any_MIS': sum(cnt.get(k, 0) for k in MIS) / tot,
            'per_MIS': {lab(k): cnt.get(k, 0) / tot for k in MIS}, 'P_MWIS': cnt.get(bb, 0) / tot, 'MWIS': lab(bb),
            'counts': {lab(k): v for k, v in cnt.items()}}

vendor = NoiseModel(p_false_pos=0.02, p_false_neg=0.10, dephasing_rate=1 / 4.5, runs=40, samples_per_run=50)
res = {}
RubyDMM = dataclasses.replace(RubyLike, name='RubyLike+DMM', dmm_objects=(DMM(bottom_detuning=-60, total_bottom_detuning=-600),))
for TT in TS:
    T = TT
    for tag, cfg in (('ideal', None), ('noisy_vendor_nominal', SimConfig.from_noise_model(vendor))):
        sim = QutipEmulator.from_sequence(seq_for(one, RubyLike), sampling_rate=0.05, config=cfg)
        res[f'T{TT}_MIS_{tag}'] = stats(sim.run().sample_final_state(N_samples=4000))
        sim = QutipEmulator.from_sequence(seq_for(one, RubyDMM, weighted=True), sampling_rate=0.05, config=cfg)
        res[f'T{TT}_MWIS_{tag}'] = stats(sim.run().sample_final_state(N_samples=4000))
        print(TT, tag, {k: round(v, 3) for k, v in res[f'T{TT}_MWIS_{tag}'].items() if k in ('P_valid', 'P_MWIS')}, flush=True)
json.dump({'register_um': coords, 'preflight': PRE, 'emulation': res}, open(f'{OUT}/emu_dist.json', 'w'), indent=1, default=float)
print('done', len(res), 'runs')
