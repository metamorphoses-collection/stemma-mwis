"""Score the Ruby Lai de l'Ombre MIS run against rules fixed before the run.
Written and self-tested before any hardware result. Do not change the rules after looking at data.

usage: python3 score_ruby_mis_lai.py lai_mis_prep.json <result.json> [...] [--order q0,q1,...]
Result shapes as in score_ruby_mis.py (device driver layout: counts under 'results', order under 'qubit_order').
"""
import json, sys, math, itertools as it
from collections import Counter

FAM = ['CG', 'DE', 'AB', 'DEF', 'EF', 'CEG', 'DF']
CONFLICT = {('DE', 'EF'), ('DE', 'CEG'), ('DE', 'DF'), ('DEF', 'CEG'), ('EF', 'CEG'), ('EF', 'DF')}
conf = lambda a, b: (a, b) in CONFLICT or (b, a) in CONFLICT
MIS = [frozenset(s) for s in ({'CG', 'DE', 'AB', 'DEF'}, {'CG', 'AB', 'DEF', 'EF'}, {'CG', 'AB', 'DEF', 'DF'}, {'CG', 'AB', 'CEG', 'DF'})]
EDITION = MIS[0]                  # unique MWIS (weight 281) = the tree of Bedier 1890 / Quentin, unrooted
NATOMS = 56
EMU = {'noiseless': {'P_valid': 1.0, 'P_any_MIS': 0.7375, 'per_MIS': {'{CG,DE,AB,DEF}': 0.28, '{CG,AB,DEF,EF}': 0.0455, '{CG,AB,DEF,DF}': 0.30925, '{CG,AB,CEG,DF}': 0.10275}},
       'noisy': {'P_valid': 0.9325, 'P_any_MIS': 0.21325, 'per_MIS': {'{CG,DE,AB,DEF}': 0.06725, '{CG,AB,DEF,EF}': 0.0345, '{CG,AB,DEF,DF}': 0.064, '{CG,AB,CEG,DF}': 0.0475}}}
PRISCUS_RATIO = 0.62 / 0.51       # Priscus rows 728-731: hardware / noisy-emulated P(any MIS)
RANDOM = 4 / 128


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
    assert len(order) == NATOMS and set(order) == set(prep['register_um']), f'qubit order must be the {NATOMS} register ids'

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
    assert all(len(k) == NATOMS and set(k) <= {'0', '1'} for k in counts), f'bitstrings must be {NATOMS} characters of 0/1'
    copies = sorted({q.split('_')[0] for q in order}, key=lambda c: int(c[1:]))
    idx = {q: i for i, q in enumerate(order)}
    per = Counter(); per_copy = {c: Counter() for c in copies}
    for bs, n in counts.items():
        for c in copies:
            on = frozenset(f for f in FAM if bs[idx[f'{c}_{f}']] == '1')
            per[on] += n; per_copy[c][on] += n
    N = sum(per.values())
    lab = lambda st: '{' + ','.join(f for f in FAM if f in st) + '}'
    valid = sum(v for s, v in per.items() if not any(conf(a, b) for a, b in it.combinations(s, 2)))
    anymis = sum(per[s] for s in MIS)
    lv, _ = wilson(valid, N); lm, _ = wilson(anymis, N)
    shares = {lab(s): per[s] / N for s in MIS}
    H1 = lv >= 0.80
    H2 = lm > RANDOM and anymis / N >= 0.10
    H3 = all(v >= 0.02 for v in shares.values())
    tv = {k: 0.5 * sum(abs(shares[m] - EMU[k]['per_MIS'][m]) for m in shares) + 0.5 * abs((1 - sum(shares.values())) - (1 - sum(EMU[k]['per_MIS'].values())))
          for k in EMU}
    out = dict(files=per_file, shots=sum(per_file.values()), copy_samples=N,
               P_valid=valid / N, P_valid_wilson_lo=lv, H1_valid_lo_ge_0_80=H1,
               P_any_MIS=anymis / N, P_any_MIS_wilson_lo=lm, H2_above_random_and_ge_0_10=H2,
               per_MIS_share=shares, H3_all_four_ge_0_02=H3,
               reported_P_edition_tree=per[EDITION] / N,
               reported_ratio_to_noisy=(anymis / N) / EMU['noisy']['P_any_MIS'], reported_priscus_ratio=PRISCUS_RATIO,
               reported_TV_MIS_shares=tv,
               reported_P_isolated_on={f: sum(v for s, v in per.items() if f in s) / N for f in ('AB', 'CG')},
               per_copy_P_any_MIS={c: sum(per_copy[c][s] for s in MIS) / max(1, sum(per_copy[c].values())) for c in copies},
               top_outcomes={lab(st): v / N for st, v in per.most_common(10)}, emulated=EMU)
    print(json.dumps(out, indent=1))
    print('\nVERDICT (frozen rules): H1', 'PASS' if H1 else 'FAIL', '| H2', 'PASS' if H2 else 'FAIL', '| H3', 'PASS' if H3 else 'FAIL',
          f'| any MIS {anymis / N:.3f} = {out["reported_ratio_to_noisy"]:.2f} x noisy emulation (Priscus {PRISCUS_RATIO:.2f})')


if __name__ == '__main__':
    main()
