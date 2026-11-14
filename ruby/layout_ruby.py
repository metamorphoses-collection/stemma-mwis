"""Re-optimise the 5-family unit-disk layout for Ruby: a larger minimum distance (relative to the blockade scale) so the
5-um spacing floor allows a tighter, stronger register. usage: python3 layout_ruby.py <layout6.json> <out.json>"""
import json, sys, numpy as np, itertools as it
from scipy.optimize import minimize
L = json.load(open(sys.argv[1])); names, W, adj = L['names'], L['W'], np.array(L['adj']); n = len(names); iu = np.triu_indices(n, 1)
def loss(x, ein=0.80, eout=1.30, dmin=0.72):
    X = x.reshape(n, 2); D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1) + 1e-12)[iu]; A = adj[iu]
    return (np.maximum(0, D - ein) ** 2 * A).sum() + (np.maximum(0, eout - D) ** 2 * (1 - A)).sum() + (np.maximum(0, dmin - D) ** 2).sum()
rng = np.random.default_rng(7); best = None
for t in range(3000):
    r = minimize(loss, rng.normal(scale=1.0, size=2 * n), method='L-BFGS-B')
    if best is None or r.fun < best.fun: best = r
    if r.fun < 1e-14: break
X = best.x.reshape(n, 2); X -= X.mean(0); D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))[iu]; A = adj[iu]
print('loss', best.fun, 'max conflict', D[A == 1].max().round(3), 'min compatible', D[A == 0].min().round(3), 'min any', D.min().round(3))
json.dump({'names': names, 'W': W, 'X': X.tolist(), 'adj': adj.tolist()}, open(sys.argv[2], 'w'))
