"""Post-hoc sensitivity check for the Lawrence sampler share (not part of the frozen analysis): whole-shot bootstrap that
keeps the seven copies of each shot together, so that correlation between copies in one shot is respected.
usage: shot_bootstrap_lawrence.py lawrence_mis_prep.json ruby_diag_job_732.json [replicates]"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import score_ruby_mis_lawrence as S
prep = json.load(open(sys.argv[1])); res = json.load(open(sys.argv[2])); B = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
order = list(prep['register_um']); assert res['qubit_order'] == order
idx = {q: i for i, q in enumerate(order)}
copies = sorted({q.split('_')[0] for q in order}, key=lambda c: int(c[1:]))
shots = []
for bs, n in res['results'].items():
    e = r = 0
    for c in copies:
        on = frozenset(f for f in S.FAM if bs[idx[f'{c}_{f}']] == '1')
        e += on == S.E; r += on == S.R
    shots += [(e, r)] * int(n)
a = np.array(shots); print('shots', len(a), 'E', a[:, 0].sum(), 'R', a[:, 1].sum(), 'share', a[:, 0].sum() / a.sum())
rng = np.random.default_rng(20261115)
ix = rng.integers(0, len(a), (B, len(a))); e = a[ix, 0].sum(1); r = a[ix, 1].sum(1); s = e / (e + r)
print('boot mean %.4f  95%% %.4f-%.4f' % (s.mean(), *np.percentile(s, [2.5, 97.5])))
