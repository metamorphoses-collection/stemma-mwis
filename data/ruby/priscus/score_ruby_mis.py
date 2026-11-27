"""Score the Ruby Priscus MIS run against rules fixed before the run.
Written before any hardware result. Do not change the rules after looking at data.

usage: python3 score_ruby_mis.py ruby_mis_prep.json <result.json> [<result2.json> ...] [--order q0,q1,...]
(Amendment 1: rows 728-731 are pooled by passing all four files; counts are summed. Rules unchanged.)

<result.json>: the Ruby result file. Accepted shapes: {"counts": {bitstring: n}}, {bitstring: n}, or
{"samples": [bitstring, ...]}. Bitstrings must follow the REGISTER QUBIT ORDER of the submitted
sequence (Pulser's convention). The default order is ruby_mis_prep.json["register_um"] key order
(c0_BMP, c0_MP, c0_EB, c0_EP, c0_EMP, c1_BMP, ...). If the submitted Register was rebuilt or relabelled,
pass the actual order with --order, or the per-copy split is wrong. ASSERT the order on the Mac before
scoring: print reg.qubit_ids from the job's sequence.
"""
import json, sys, math, itertools as it
from collections import Counter

FAM = ['BMP', 'MP', 'EB', 'EP', 'EMP']
CONFLICT = {('BMP', 'EB'), ('BMP', 'EP'), ('BMP', 'EMP'), ('MP', 'EP'), ('EB', 'EP'), ('EB', 'EMP')}
conf = lambda a, b: (a, b) in CONFLICT or (b, a) in CONFLICT
MIS = [frozenset(s) for s in ({'BMP', 'MP'}, {'MP', 'EB'}, {'MP', 'EMP'}, {'EP', 'EMP'})]
DEBOOR = frozenset({'BMP', 'MP'})
EMU_NOISY = {'{BMP,MP}': 0.124, '{MP,EB}': 0.108, '{MP,EMP}': 0.166, '{EP,EMP}': 0.114}   # frozen prediction
RANDOM_BASELINE = 4 / 32

def wilson(k, n, z=1.96):
    if n == 0: return (0, 0)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h

args = [a for a in sys.argv[1:]]
order_arg = None
if '--order' in args:
    i = args.index('--order'); order_arg = args[i + 1]; del args[i:i + 2]
prep = json.load(open(args[0])); files = args[1:]
assert files, 'give at least one result file'
order = list(prep['register_um'])
if order_arg: order = order_arg.split(',')
assert len(order) == 45 and set(order) == set(prep['register_um']), 'qubit order must be the 45 register ids'
def load_counts(path):
    res = json.load(open(path))
    c = res.get('counts') if isinstance(res, dict) and 'counts' in res else None
    if c is None and isinstance(res, dict) and 'results' in res and isinstance(res['results'], dict):
        c = res['results']  # device driver layout (ruby_diag_job_NNN.json): counts under 'results', atom order under 'qubit_order'
        assert res.get('qubit_order') == order, 'file qubit_order differs from the scoring order'
    if c is None and isinstance(res, dict) and 'samples' in res: c = Counter(res['samples'])
    if c is None: c = res
    return {str(k): int(v) for k, v in c.items()}
counts = Counter()
per_file_shots = {}
for f in files:
    c = load_counts(f); per_file_shots[f] = sum(c.values()); counts.update(c)
counts = dict(counts)
assert all(len(k) == 45 and set(k) <= {'0', '1'} for k in counts), 'bitstrings must be 45 characters of 0/1'
shots = sum(counts.values())

copies = sorted({q.split('_')[0] for q in order}, key=lambda c: int(c[1:]))
idx = {q: i for i, q in enumerate(order)}
per = Counter(); per_copy = {c: Counter() for c in copies}
for bs, nshot in counts.items():
    for c in copies:
        on = frozenset(f for f in FAM if bs[idx[f'{c}_{f}']] == '1')
        per[on] += nshot; per_copy[c][on] += nshot
N = sum(per.values())
valid = sum(v for s, v in per.items() if not any(conf(a, b) for a, b in it.combinations(s, 2)))
anymis = sum(per[s] for s in MIS)
lab = lambda s: '{' + ','.join(f for f in FAM if f in s) + '}'
shares = {lab(s): per[s] / N for s in MIS}
tv = 0.5 * sum(abs(shares[k] - EMU_NOISY[k]) for k in EMU_NOISY) + 0.5 * max(0.0, 1 - sum(shares.values()) - (1 - sum(EMU_NOISY.values())))

lv, _ = wilson(valid, N); lm, _ = wilson(anymis, N)
H1 = lv >= 0.80
H2 = (lm > RANDOM_BASELINE) and (anymis / N >= 0.25)
H3 = all(per[s] / N >= 0.02 for s in MIS)
out = {'files': per_file_shots, 'shots': shots, 'copy_samples': N,
       'P_valid': valid / N, 'P_valid_wilson_lo': lv, 'H1_valid_ge_0.80': H1,
       'P_any_MIS': anymis / N, 'P_any_MIS_wilson_lo': lm, 'H2_above_random_and_ge_0.25': H2,
       'per_MIS_share': shares, 'H3_all_four_ge_0.02': H3,
       'P_deBoor_tree': per[DEBOOR] / N,
       'reported_TV_vs_noisy_emulation': tv,
       'per_copy_P_any_MIS': {c: sum(per_copy[c][s] for s in MIS) / max(1, sum(per_copy[c].values())) for c in copies},
       'top_outcomes': {lab(s): v / N for s, v in per.most_common(8)}}
print(json.dumps(out, indent=1))
print('\nVERDICT (frozen rules):', 'H1', 'PASS' if H1 else 'FAIL', '| H2', 'PASS' if H2 else 'FAIL', '| H3', 'PASS' if H3 else 'FAIL')
