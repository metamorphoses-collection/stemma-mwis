"""Find a 2D layout for a conflict graph at Ruby-like geometry: conflicts <= c, compatible >= 1.47c (residual tail),
all pairs >= dmin (layout units; the prep script scales by 7.2 um/unit, so c=0.90 -> 6.5 um < R_b 6.73 um,
1.27c -> 8.2 um >= 1.25 R_b, dmin 0.70 -> 5.0 um). Restarts until zero violation. Not a proof of impossibility.
Usage: layout_hw.py graph.json out_layout.json [restarts]"""
import sys, json, numpy as np
from scipy.optimize import minimize
g = json.load(open(sys.argv[1])); n = g['n']; A = np.array(g['adj']); iu = np.triu_indices(n, 1); a = A[iu]
C, RAT, DMIN = 0.90, 1.47, 0.72
def loss(x):
    X = x.reshape(n, 2); D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1) + 1e-12)[iu]
    return ((np.maximum(0, D - C) ** 2) * a).sum() + ((np.maximum(0, RAT * C - D) ** 2) * (1 - a)).sum() \
        + (np.maximum(0, DMIN - D) ** 2).sum()
rng = np.random.default_rng(1); best = None
for t in range(int(sys.argv[3]) if len(sys.argv) > 3 else 2000):
    r = minimize(loss, rng.normal(scale=np.sqrt(n) / 2, size=2 * n), method='L-BFGS-B', options={'maxiter': 5000})
    if best is None or r.fun < best.fun: best = r
    if r.fun < 1e-14: break
X = best.x.reshape(n, 2); D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))[iu]
viol = int(((D > C + 1e-9) & (a == 1)).sum() + ((D < RAT * C - 1e-9) & (a == 0)).sum() + (D < DMIN).sum())
print(f'n={n} edges={int(a.sum())} restarts={t+1} loss={best.fun:.3g} violations={viol}',
      f'max conflict {D[a==1].max():.3f} min compatible {D[a==0].min():.3f} min pair {D.min():.3f}')
if viol == 0:
    json.dump(dict(names=g['names'], W=g['W'], X=X.tolist(), adj=g['adj']), open(sys.argv[2], 'w'))
    print('layout written')
