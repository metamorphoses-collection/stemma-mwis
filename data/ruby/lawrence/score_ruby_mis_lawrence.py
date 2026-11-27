"""Score the Ruby Lawrence-core MIS run against rules fixed before the run.
Written and self-tested before any hardware result. Do not change the rules after looking at data.

usage: python3 score_ruby_mis_lawrence.py lawrence_mis_prep.json <result.json> [...] [--order q0,q1,...]
Result shapes as in score_ruby_mis.py (device driver layout: counts under 'results', order under 'qubit_order').
"""
import json, sys, math, itertools as it
from collections import Counter

FAM = ['B+F+Mi+P+R', 'B+R', 'P+R', 'B+F+P+R', 'B+F+R']
CONFLICT = {('P+R', 'B+R'), ('P+R', 'B+F+R')}
conf = lambda a, b: (a, b) in CONFLICT or (b, a) in CONFLICT
E = frozenset({'B+F+Mi+P+R', 'B+R', 'B+F+P+R', 'B+F+R'})      # edition's history = unique MIS = MWIS
R = frozenset({'B+F+Mi+P+R', 'P+R', 'B+F+P+R'})               # rival history
BOOT = 0.795; BAND = (0.745, 0.845)
EMU = {'noiseless': {'P_valid': 0.99975, 'P_E': 0.685, 'P_R': 0.032, 'share': 0.956},
       'noisy': {'P_valid': 0.95175, 'P_E': 0.140, 'P_R': 0.071, 'share': 0.664}}
RANDOM = 1 / 32


def wilson(k, n, z=1.96):
    if n == 0: return (0.0, 0.0)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def main():
    args = sys.argv[1:]; order_arg = None
    if '--order' in args:
        i = args.index('--order'); order_arg = args[i + 1]; del args[i:i + 2]
    prep = json.load(open(args[0])); files = args[1:]
    assert files, 'give at least one result file'
    order = order_arg.split(',') if order_arg else list(prep['register_um'])
    assert len(order) == 35 and set(order) == set(prep['register_um']), 'qubit order must be the 35 register ids'

    def load(path):
        res = json.load(open(path))
        if isinstance(res, dict) and 'results' in res and isinstance(res['results'], dict):
            assert res.get('qubit_order') == order, 'file qubit_order differs from the scoring order'
            c = res['results']
        elif isinstance(res, dict) and 'counts' in res: c = res['counts']
        elif isinstance(res, dict) and 'samples' in res: c = Counter(res['samples'])
        else: c = res
        return {str(k): int(v) for k, v in c.items()}

    counts = Counter(); per_file = {}
    for f in files:
        c = load(f); per_file[f] = sum(c.values()); counts.update(c)
    assert all(len(k) == 35 and set(k) <= {'0', '1'} for k in counts), 'bitstrings must be 35 characters of 0/1'
    copies = sorted({q.split('_')[0] for q in order}, key=lambda c: int(c[1:]))
    idx = {q: i for i, q in enumerate(order)}
    per = Counter(); per_copy = {c: Counter() for c in copies}
    for bs, n in counts.items():
        for c in copies:
            on = frozenset(f for f in FAM if bs[idx[f'{c}_{f}']] == '1')
            per[on] += n; per_copy[c][on] += n
    N = sum(per.values())
    valid = sum(v for s, v in per.items() if not any(conf(a, b) for a, b in it.combinations(s, 2)))
    nE, nR = per[E], per[R]; pE, pR = nE / N, nR / N
    lv, _ = wilson(valid, N); lE, _ = wilson(nE, N)
    sd = math.sqrt(max(pE + pR - (pE - pR) ** 2, 0) / N); dlo = (pE - pR) - 1.96 * sd
    H1 = lv >= 0.80
    H2 = lE > RANDOM and pE >= 0.07
    H3 = dlo > 0
    s = nE / (nE + nR) if nE + nR else float('nan'); slo, shi = wilson(nE, nE + nR)
    S = 'match' if (BAND[0] <= slo and shi <= BAND[1]) else ('mismatch' if (shi < BAND[0] or slo > BAND[1]) else 'undetermined')
    lab = lambda st: '{' + ','.join(f for f in FAM if f in st) + '}'
    out = dict(files=per_file, shots=sum(per_file.values()), copy_samples=N,
               P_valid=valid / N, P_valid_wilson_lo=lv, H1_valid_lo_ge_0_80=H1,
               P_E=pE, P_E_wilson_lo=lE, H2_E_above_random_and_ge_0_07=H2,
               P_R=pR, diff_E_minus_R=pE - pR, diff_lo=dlo, H3_E_more_than_R=H3,
               share=s, share_ci=[slo, shi], S_sampler=S, bootstrap=BOOT, band=BAND, emulated=EMU,
               per_copy_P_E={c: per_copy[c][E] / max(1, sum(per_copy[c].values())) for c in copies},
               top_outcomes={lab(st): v / N for st, v in per.most_common(10)})
    print(json.dumps(out, indent=1))
    print('\nVERDICT (frozen rules): H1', 'PASS' if H1 else 'FAIL', '| H2', 'PASS' if H2 else 'FAIL', '| H3', 'PASS' if H3 else 'FAIL',
          f'| S {S} (share {s:.3f}, CI {slo:.3f}-{shi:.3f}; bootstrap {BOOT})')


if __name__ == '__main__':
    main()
