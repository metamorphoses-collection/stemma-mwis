r"""How much a localised margin certificate could gain on the Lai (Discussion of the technical report).

Theorem gammaconf pays for all N = 56 splits on 7 witnesses through Clopper-Pearson bounds at level delta/(2N).
A localised certificate would pay for fewer families, which can only narrow the intervals; the most favourable
case is N_paid = 1. For N_paid in {56, 8, 4, 1} this script finds the smallest multiple k of the Lai's length
(same proportions, step 0.25) at which gamma_lo of the estimated stemma becomes positive.
Usage: localised_bound.py OUT.json"""
import sys, json
src = open(__file__.replace('localised_bound.py', 'gamma_confidence.py')).read().split('\n')
ns = {}
exec('\n'.join(src[:135]), ns)          # imports, helpers, Solver, gamma_of, cp_bounds
exec('\n'.join(src[175:193]), ns)       # Lai units, counts, estimated stemma T
cnt, L, T, sol7, cp_bounds, gamma_of = (ns[k] for k in ('cnt', 'L', 'T', 'sol7', 'cp_bounds', 'gamma_of'))
D = 0.05
def glo(k, n_paid):
    c2 = {v: int(round(k * x)) for v, x in cnt.items()}; L2 = int(round(k * L)); d = D / (2 * n_paid)
    b = {v: cp_bounds(c2.get(v, 0), L2, d) for v in sol7.V}
    return gamma_of(sol7, T, {v: b[v][0] for v in sol7.V}, {v: b[v][1] for v in sol7.V})
out = {}
for n_paid in (56, 8, 4, 1):
    ks = [k / 4 for k in range(4, 25)]
    first = next((k for k in ks if glo(k, n_paid) > 0), None)
    out[n_paid] = dict(gamma_lo_at_1=glo(1, n_paid), first_multiple_certified=first)
    print(n_paid, out[n_paid])
json.dump(out, open(sys.argv[1], 'w'), indent=1)
assert out[1]['gamma_lo_at_1'] < 0 and out[1]['first_multiple_certified'] > 1
