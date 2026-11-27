"""Lai de l'Ombre: unit-disk layout as layout_ruby.py (same margins: conflict <= 0.80, compatible >= 1.30, any >= 0.72
layout units of 7.2 um), plus a compactness term so more copies tile inside 50 um. Margins re-checked after.
usage: python3 lai_layout_compact.py <layout_in.json> <out.json>"""
import json, sys, numpy as np
from scipy.optimize import minimize
L = json.load(open(sys.argv[1])); names, W, adj = L['names'], L['W'], np.array(L['adj']); n = len(names); iu = np.triu_indices(n, 1)
def hinge(x, ein=0.79, eout=1.31, dmin=0.73):
    X = x.reshape(n, 2); D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1) + 1e-12)[iu]; A = adj[iu]
    return (np.maximum(0, D - ein) ** 2 * A).sum() + (np.maximum(0, eout - D) ** 2 * (1 - A)).sum() + (np.maximum(0, dmin - D) ** 2).sum()
def rmax(x, p=12):
    X = x.reshape(n, 2); X = X - X.mean(0); return ((((X ** 2).sum(1)) ** (p / 2)).sum()) ** (1 / p)
rng = np.random.default_rng(7); best = None
for t in range(400):
    r = minimize(lambda x: 1e4 * hinge(x) + rmax(x), rng.normal(scale=1.0, size=2 * n), method='L-BFGS-B')
    if hinge(r.x) < 1e-10 and (best is None or rmax(r.x) < rmax(best)): best = r.x
X = best.reshape(n, 2); X -= X.mean(0); D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))[iu]; A = adj[iu]
assert D[A == 1].max() <= 0.80 and D[A == 0].min() >= 1.30 and D.min() >= 0.72, 'margins violated'
print('max conflict', D[A == 1].max().round(3), 'min compatible', D[A == 0].min().round(3), 'min any', D.min().round(3), 'radius um', (np.sqrt((X ** 2).sum(1)).max() * 7.2).round(2))
json.dump({'names': names, 'W': W, 'X': X.tolist(), 'adj': adj.tolist()}, open(sys.argv[2], 'w'))
