"""Task 3: weighted and unweighted emulations, single copy, the preps' drive and noise model, T and the weight floor varied.
usage: python sampler_test_emulate.py <instance> <layout> [T list]
  instance: priscus | lai | lawrence ; layout: register (frozen c0 copy) | opt<R> (sampler_test_layout_lawrence.json,
  footprint radius R) ; T list default 1000,2000,4000,6000 (ns).
Settings per T: mode in {unweighted, w0.3, w0.1, w0.0} x noise in {ideal, noisy}. a = 0.0 is an extra (pure
proportional weights, no per-family floor), not in the brief.
Per setting: 4000 samples (np.random seeded per setting) via sample_final_state (readout errors applied by Pulser in
the noisy case), PLUS the exact outcome distribution: |amplitude|^2 (ideal) or diag(rho) (noisy) pushed through the
same independent readout channel (0->1 with 0.02, 1->0 with 0.10). The exact distribution is checked against the
samples (total-variation distance reported).
Output: sampler_test_20261123/emu_<instance>_<layout>.json (rewritten after every setting)."""
import sys, os, json, time, dataclasses, warnings, itertools as it
import numpy as np
from pulser import Register, Sequence, Pulse, NoiseModel
from pulser.waveforms import InterpolatedWaveform, RampWaveform
from pulser.devices import AnalogDevice
from pulser.channels.dmm import DMM
from pulser_simulation import QutipEmulator, SimConfig
from sampler_test_common import INST, OUT, OMEGA, DELTA0, DELTAF, C6, register_positions, weff, geometry, ground_check, lab

warnings.filterwarnings('ignore')
name, layout = sys.argv[1], sys.argv[2]
TS = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else '1000,2000,4000,6000').split(',')]
g = INST[name]; F = g['fam']; n = len(F)
assert abs(AnalogDevice.interaction_coeff - C6) < 1
if layout == 'register':
    P = register_positions(name)
else:
    L = json.load(open(os.path.join(OUT, f'layout_{name}.json')))['optimised'][layout[3:]]
    P = np.array([L['positions_um'][f] for f in F], float)
reg = Register({f: tuple(map(float, P[i])) for i, f in enumerate(F)})
Dev = dataclasses.replace(AnalogDevice.to_virtual(), name='RubyLike+DMM', max_radial_distance=50, max_atom_num=100,
                          min_atom_distance=5, dmm_objects=(DMM(bottom_detuning=-60, total_bottom_detuning=-600),))
OFF_DEVICE = max(TS) > 6000
if OFF_DEVICE:   # off-device adiabaticity check only: Ruby's 6000 ns cap lifted
    Dev = dataclasses.replace(Dev, max_sequence_duration=None)
FP, FN = 0.02, 0.10
vendor = NoiseModel(p_false_pos=FP, p_false_neg=FN, dephasing_rate=1 / 4.5, runs=40, samples_per_run=50)


def seq(T, a):
    s = Sequence(reg, Dev); s.declare_channel('ryd', 'rydberg_global')
    if a is not None:
        we = weff(g, a)
        s.config_detuning_map(reg.define_detuning_map({q: 1 - we[q] for q in reg.qubit_ids}, 'dmap'), 'dmm_0')
    s.add(Pulse(InterpolatedWaveform(T, [1e-9, OMEGA, OMEGA, 1e-9]),
                InterpolatedWaveform(T, [DELTA0, DELTA0 * 0.5, DELTAF * 0.6, DELTAF]), 0), 'ryd')
    if a is not None:
        s.add_dmm_detuning(RampWaveform(T, 0, -DELTAF), 'dmm_0')
    return s


BITS = [''.join(b) for b in it.product('01', repeat=n)]


def exact_dist(state, noisy):
    """Outcome distribution over bitstrings ('1' = Rydberg). Pulser ground-rydberg basis: index 0 = |r>, 1 = |g>."""
    M = state.full()
    p = np.real(np.diag(M)) if M.shape[0] == M.shape[1] else np.abs(M[:, 0]) ** 2
    pre = {}
    for k, v in enumerate(p):
        b = ''.join('1' if ch == '0' else '0' for ch in format(k, f'0{n}b')); pre[b] = float(v)
    if not noisy: return pre
    post = {b: 0.0 for b in BITS}
    for b, v in pre.items():
        if v < 1e-15: continue
        for o in BITS:
            pr = 1.0
            for x, y in zip(b, o):
                pr *= (1 - FN if y == '1' else FN) if x == '1' else (FP if y == '1' else 1 - FP)
            post[o] += v * pr
    return post


SUFFIX = sys.argv[4] if len(sys.argv) > 4 else ''   # e.g. '_longT' for the off-device T scan
outp = os.path.join(OUT, f'emu_{name}_{layout}{SUFFIX}.json')
res = dict(instance=name, layout=layout, positions_um={f: P[i].tolist() for i, f in enumerate(F)}, geometry=geometry(g, P),
           ground_checks={f'a{a}': ground_check(g, P, a) for a in (0.3, 0.1, 0.0)},
           drive=dict(Omega=OMEGA, delta0=DELTA0, deltaf=DELTAF, shape='Omega 4 knots 0,O,O,0; delta 4 knots d0,d0/2,0.6df,df; '
                      'DMM map 1-w_eff, ramp 0 -> -df'), noise=dict(p_false_pos=FP, p_false_neg=FN, dephasing_rate=1 / 4.5),
           samples=4000, sampling_rate=0.05, off_device_duration=OFF_DEVICE, engine='pulser 1.9.1 QutipEmulator, qutip 5.3.1', runs={})
seed = 0
for T in TS:
    for mode, a in (('unweighted', None), ('w0.3', 0.3), ('w0.1', 0.1), ('w0.0', 0.0)):
        for noise in ('ideal', 'noisy'):
            seed += 1; t0 = time.time()
            sim = QutipEmulator.from_sequence(seq(T, a), sampling_rate=0.05,
                                              config=SimConfig.from_noise_model(vendor) if noise == 'noisy' else None)
            r = sim.run()
            np.random.seed(1000 * TS.index(T) + seed)
            cnt = r.sample_final_state(N_samples=4000)
            ex = exact_dist(r.get_final_state(), noise == 'noisy')
            emp = {b: cnt.get(b, 0) / 4000 for b in BITS}
            tv = 0.5 * sum(abs(emp[b] - ex[b]) for b in BITS)
            key = f'T{T}_{mode}_{noise}'
            res['runs'][key] = dict(T=T, mode=mode, a=a, noise=noise, seconds=time.time() - t0, tv_samples_vs_exact=tv,
                                    counts={lab(g, {F[i] for i in range(n) if b[i] == '1'}): int(v) for b, v in cnt.items()},
                                    exact={lab(g, {F[i] for i in range(n) if b[i] == '1'}): v for b, v in ex.items() if v > 1e-9})
            json.dump(res, open(outp, 'w'), indent=1, default=float)
            print(key, f'{time.time() - t0:.1f}s', f'TV={tv:.3f}', sorted(cnt.items(), key=lambda t: -t[1])[:3], flush=True)
print('done', outp)
