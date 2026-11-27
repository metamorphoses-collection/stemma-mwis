"""Q3/Q4/Q5/Q8 support. Reads the frozen emulation outputs (no new emulation is run) and
(1) tabulates emulated predictions with their 4000-sample binomial SE (conflict-free, maximum, Lawrence E/(E+R));
(2) Lawrence E/(E+R) across the 16 emulated settings (data/lawrence_emu_dist.json) and the size-bias reading:
    under P(S) ~ lambda^|S| over independent sets, E/(E+R) = lambda/(1+lambda) (E has 4 families, R 3);
    and the bootstrap condition (E beats R iff w(B+R)+w(B+F+R) > w(P+R));
(3) Q8: classical (Omega=0) Rydberg energy at the end of the sweep WITH all van der Waals tails, exhaustively over the
    2^n configurations of one copy, unweighted (global delta_f) and weighted (delta_i = delta_f * w_eff_i, the DMM map of
    ruby_mis_prep.py); does its minimum coincide with the abstract-graph MIS / MWIS? Plus the summed inter-copy tail.
Output: checks_20261123/emulation_and_tails.json"""
import json, os, math, itertools as it
import numpy as np
from check_common import HERE, GRAPHS, LAW_E, LAW_R, independent_sets, maximum_sets, lab, dump

C6 = 865723.02; DELTAF = 2 * math.pi * 0.9
STEM = os.path.dirname(os.path.dirname(HERE))
out = {'emulation_settings': dict(
    engine='pulser_simulation.QutipEmulator, single copy of the instance, sampling_rate=0.05',
    drive='Omega 0->2pi*1.5->2pi*1.5->0 rad/us; delta -2pi*4 -> -2pi*2 -> 2pi*0.54 -> 2pi*0.9 rad/us (interpolated); T 6000 ns',
    samples='sample_final_state(N_samples=4000) per configuration',
    noiseless='pure-state Schroedinger evolution, no noise',
    noisy='NoiseModel(p_false_pos=0.02, p_false_neg=0.10, dephasing_rate=1/4.5 us^-1): readout errors + Lindblad dephasing, '
          'solved as a density matrix (master equation). runs=40, samples_per_run=50 were passed but Pulser warns they are NOT '
          'used by the active noise types (SPAM, dephasing) [data/lai_ruby_dist/log.txt] -> no trajectories; 4000 samples from one density matrix.',
    weighted='DMM detuning map (1 - w_eff), ramp 0 -> -delta_f over T; w_eff = 0.3 + 0.7 w/w_max; device = AnalogDevice virtual + DMM(bottom -60, total -600)')}

se4k = lambda p: math.sqrt(p * (1 - p) / 4000)
preds = {}
pr = json.load(open(os.path.join(HERE, 'ruby_mis_prep.json')))['emulation']
lai = json.load(open(os.path.join(STEM, '..', '..', '..', 'code', 'stemma-mwis', 'data', 'ruby', 'lai', 'lai_emu_dist_T6000.json')))['emulation'] \
    if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lai', 'lai_emu_dist_T6000.json')) else None
lai = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lai', 'lai_emu_dist_T6000.json')))['emulation']
laiprep = json.load(open(os.path.join(HERE, 'lai_mis_prep.json')))['emulation']
law = json.load(open(os.path.join(STEM, 'data', 'lawrence_emu_dist_T6000.json')))['emulation']
lawscan = json.load(open(os.path.join(STEM, 'data', 'lawrence_emu_dist.json')))['emulation']
for tag, src in (('priscus', {'ideal': pr['MIS_ideal'], 'noisy': pr['MIS_noisy_vendor_nominal'],
                              'MWIS_ideal': pr['MWIS_ideal'], 'MWIS_noisy': pr['MWIS_noisy_vendor_nominal']}),
                 ('lai (scorer, lai_emu_dist_T6000.json)', {'ideal': lai['T6000_MIS_ideal'], 'noisy': lai['T6000_MIS_noisy_vendor_nominal']}),
                 ('lai (lai_mis_prep.json, independent rerun)', {'ideal': laiprep['MIS_ideal'], 'noisy': laiprep['MIS_noisy_vendor_nominal']}),
                 ('lawrence (lawrence_emu_dist_T6000.json)', {'ideal': law['T6000_MIS_ideal'], 'noisy': law['T6000_MIS_noisy_vendor_nominal']})):
    preds[tag] = {k: dict(P_valid=v['P_valid'], P_valid_se=se4k(v['P_valid']), P_any_MIS=v['P_any_MIS'], P_any_MIS_se=se4k(v['P_any_MIS']),
                          P_MWIS=v.get('P_MWIS')) for k, v in src.items()}
out['predictions'] = preds

# Lawrence E/(E+R) in emulation
lE, lR = '{B+F+Mi+P+R,B+R,B+F+P+R,B+F+R}', '{B+F+Mi+P+R,P+R,B+F+P+R}'
def share(c):
    e, r = c.get(lE, 0), c.get(lR, 0); s = e / (e + r); return dict(E=e, R=r, share=s, se=math.sqrt(s * (1 - s) / (e + r)))
out['lawrence_share_T6000_frozen'] = {k: share(v['counts']) for k, v in law.items()}
out['lawrence_share_scan16'] = {k: share(v['counts']) for k, v in lawscan.items()}
hw = 0.7879
lam = hw / (1 - hw)
out['lawrence_size_bias'] = dict(
    E_size=4, R_size=3, E_weight=21, R_weight=18,
    bootstrap_condition='E is the MWIS iff w(B+R) + w(B+F+R) > w(P+R) (shared families cancel); full data 4+3 > 4',
    bootstrap_value=0.7955, bootstrap_file='data/lawrence_native_graph_s3_bootstrap.json (B=20000, seed 1), code/native/bootstrap_histories.py',
    hardware_share=hw, fugacity_matching_hardware=lam,
    note='device is unweighted: it can only see |E|=4 vs |R|=3. Greedy random-order gives 2/3; uniform over IS gives 1/2; '
         'a size-biased sampler P~lambda^|S| gives lambda/(1+lambda), any value in (1/2,1).',
    readout_false_negative_effect='with p_fn=0.10 per atom, E (4 atoms) survives readout with 0.9^4, R with 0.9^3: '
                                  'pre-readout share = %.3f' % ((hw / 0.9 ** 4) / (hw / 0.9 ** 4 + (1 - hw) / 0.9 ** 3)))

# Q8: tails
tails = {}
for name, g in GRAPHS.items():
    prep = json.load(open(os.path.join(HERE, g['prep'])))
    reg = prep['register_um']; weff = prep['preflight']['WEFF']; F = g['fam']; n = len(F)
    P = np.array([reg[f'c0_{f}'] for f in F]); D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf)
    V = C6 / D ** 6
    allP = np.array(list(reg.values())); ids = list(reg)
    c0 = [i for i, q in enumerate(ids) if q.startswith('c0_')]; oth = [i for i, q in enumerate(ids) if not q.startswith('c0_')]
    Dx = np.sqrt(((allP[c0][:, None] - allP[oth][None]) ** 2).sum(-1))
    inter_sum = (C6 / Dx ** 6).sum(1).max()
    ind = set(independent_sets(g)); mis = set(maximum_sets(g))
    mwis = max(ind, key=lambda s: sum(g['W'][x] for x in s)); mwis_eff = max(ind, key=lambda s: sum(weff[x] for x in s))
    rec = {'max_tail_between_compatible_um': None, 'sum_intercopy_tail_max_rad_us': float(inter_sum)}
    for mode, delta in (('unweighted', np.full(n, DELTAF)), ('weighted', np.array([DELTAF * weff[f] for f in F]))):
        Es = []
        for bits in it.product([0, 1], repeat=n):
            z = np.array(bits); E = 0.5 * z @ V @ z - delta @ z
            Es.append((E, frozenset(F[i] for i in range(n) if bits[i])))
        Es.sort(key=lambda t: t[0])
        gmin = Es[0][0]; ground = [s for e, s in Es if abs(e - gmin) < 1e-9]
        first_ex = next(e for e, s in Es if abs(e - gmin) > 1e-9)
        # lowest-energy configuration that is NOT independent
        best_bad = next((e, s) for e, s in Es if s not in ind)
        rec[mode] = dict(ground=[lab(g, s) for s in ground], ground_energy=gmin, gap_to_next=first_ex - gmin,
                         ground_equals=('MIS set' if set(ground) <= mis else 'other') if mode == 'unweighted' else
                                       ('MWIS' if ground == [mwis] else 'NOT the MWIS'),
                         n_degenerate_MIS_split_by_tails=(len(ground) < len(mis)) if mode == 'unweighted' else None,
                         MIS_energies={lab(g, s): float(0.5 * np.array([x in s for x in F]) @ V @ np.array([x in s for x in F]) - delta @ np.array([x in s for x in F])) for s in mis} if mode == 'unweighted' else None,
                         lowest_non_independent=[lab(g, best_bad[1]), best_bad[0]])
    rec['abstract_MWIS_true_w'] = lab(g, mwis); rec['abstract_MWIS_compressed_w'] = lab(g, mwis_eff)
    tails[name] = rec
out['tails_check'] = tails
print(json.dumps({k: out[k] for k in ('lawrence_share_T6000_frozen', 'lawrence_size_bias', 'tails_check')}, indent=1, default=float))
print('scan16 shares:', {k: round(v['share'], 3) for k, v in out['lawrence_share_scan16'].items()})
print('->', dump(out, 'emulation_and_tails.json'))
