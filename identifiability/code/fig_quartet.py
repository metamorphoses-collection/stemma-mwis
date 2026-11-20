"""Figure: quartet recovery vs (long pendant rate p, internal rate r), with the exact boundary b = a.
Usage: fig_quartet.py data/partA_quartet.json figures/fig_quartet.pdf"""
import sys, json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = json.load(open(sys.argv[1])); g = d['grid']
P = sorted({x['p'] for x in g}); Rr = sorted({x['r'] for x in g})
Z = np.array([[next(x['recovery'] for x in g if x['p'] == p and x['r'] == r) for p in P] for r in Rr])
fig, ax = plt.subplots(figsize=(4.6, 3.4))
im = ax.pcolormesh(P, Rr, Z, shading='nearest', cmap='Greys_r', vmin=0, vmax=1)
pp = np.linspace(min(P), max(P), 400); ps = d['ps']
# exact boundary m_true = m_false (Proposition, p1 = p3 = p, p2 = p4 = ps, r1 = r2 = r), solved for r by bisection
def gap(p, r):
    rb, pb, sb = 1 - r, 1 - p, 1 - ps
    mt = 2 * rb * pb * sb * (r + rb * p * ps)
    mf = rb * rb * (p * p * sb * sb + ps * ps * pb * pb)
    return mt - mf
rb = []
for p in pp:
    lo, hi = 1e-6, 0.5
    for _ in range(60):
        mid = (lo + hi) / 2; lo, hi = (lo, mid) if gap(p, mid) > 0 else (mid, hi)
    rb.append(hi)
ax.plot(pp, rb, color='tab:red', lw=1.5, label='exact boundary $m_{12|34}=m_{13|24}$')
ax.set_yscale('log'); ax.set_xlabel('rate of corruption of witnesses 1 and 3, $p$')
ax.set_ylabel('internal innovation rate, $r$'); ax.set_ylim(min(Rr), max(Rr))
ax.legend(loc='upper left', fontsize=7, frameon=False)
cb = fig.colorbar(im, ax=ax); cb.set_label(f"recovery, $L={d['L']}$, {d['reps']} replicates", fontsize=8)
fig.tight_layout(); fig.savefig(sys.argv[2])
