"""Check C24: the numbers stated in Section birthdeath for the pruning rule, asserted on the saved outputs to the
printed precision. Usage: bedier_calibrated_check.py EQUAL.json UNEVEN.json; exits non-zero on failure."""
import sys, json
fails = []
def check(ok, msg):
    print(('PASS ' if ok else 'FAIL ') + msg)
    if not ok: fails.append(msg)
E, U = (json.load(open(f)) for f in sys.argv[1:3])
for R, tag in ((E, 'equal'), (U, 'uneven')):
    check(sorted(k for k in R if k.startswith('alpha')) == ['alpha=0.3', 'alpha=0.6', 'alpha=0.9'], f'{tag}: three alphas')
    check(all(R[a]['true']['traditions'] == 300 for a in R if a.startswith('alpha')), f'{tag}: 300 traditions each')
fr = lambda r, a: r['recovery'][a]['false_added'] / r['true']['traditions']
tr = lambda r, a: r['recovery'][a]['true_found'] / r['recovery'][a]['true_total']
inr = lambda x, lo, hi: lo - 0.005 <= x <= hi + 0.005
FC = 'filtered calibrated'
ev = [E[a] for a in ('alpha=0.3', 'alpha=0.6', 'alpha=0.9')]
check(all(inr(fr(r, 'raw'), 0.97, 1.38) for r in ev), 'equal: raw false 0.97-1.38')
check(all(inr(fr(r, FC), 0.10, 0.19) for r in ev), 'equal: rule false 0.10-0.19')
check(all(inr(fr(r, 'k>=3'), 0.44, 0.75) for r in ev), 'equal: k>=3 false 0.44-0.75')
check(all(inr(tr(r, FC), 0.95, 1.00) for r in ev), 'equal: rule keeps 95-100%')
check(all(inr(tr(r, 'k>=3'), 0.97, 1.00) for r in ev), 'equal: k>=3 keeps 97-100%')
check(all(inr(r[FC]['bifid_root'], 0.61, 0.78) and inr(r[FC]['i3'], 0.83, 0.89) for r in ev), 'equal: rule root 0.61-0.78, i3 0.83-0.89')
uv = {a: U[a] for a in ('alpha=0.3', 'alpha=0.6', 'alpha=0.9')}
for a, (raw, rule) in {'alpha=0.3': (1.08, 0.20), 'alpha=0.6': (1.18, 0.35), 'alpha=0.9': (1.63, 0.70)}.items():
    check(abs(fr(uv[a], 'raw') - raw) < 0.005 and abs(fr(uv[a], FC) - rule) < 0.005, f'uneven {a}: raw {raw} -> rule {rule}')
check(all(inr(tr(r, FC), 0.88, 0.97) for r in uv.values()), 'uneven: rule keeps 88-97%')
fb = [R[a]['fallback_filtered'] / 300 for R in (E, U) for a in R if a.startswith('alpha')]
check(all(inr(x, 0.15, 0.19) for x in fb), 'fallback in 15-19% of traditions')
check(all(inr(U[a]['units_kept_after_filter_median'], 0.91, 0.99) or U[a]['units_kept_after_filter_median'] > 0.99 for a in uv), 'places kept 91-99%')
check(all(inr(tr(U[a], 'filtered oracle'), 0.86, 0.93) for a in uv), 'uneven: latent-parameter rule keeps 86-93%')
print('ALL PASS' if not fails else f'{len(fails)} FAILED'); sys.exit(1 if fails else 0)
