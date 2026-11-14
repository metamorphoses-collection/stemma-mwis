"""Lai de l'Ombre: retile the single-copy register of ruby_mis_prep.py on a hexagonal grid (the ring tiling fits only 4
copies of this 7-atom graph). Keeps the copy geometry, drive and preflight of the prep file; re-checks every Ruby
constraint on the tiled register. PREPARATION ONLY, submits nothing.
usage: python3 lai_tile.py <prep_in.json> <prep_out.json>"""
import json, sys, itertools as it, numpy as np
P = json.load(open(sys.argv[1])); reg = P['register_um']; C6 = 865723.02
fam = [q.split('_', 1)[1] for q in reg if q.startswith('c0_')]
one = np.array([reg[f'c0_{f}'] for f in fam]); one -= one.mean(0)
GAP = 14.0; best = None
for rot in np.linspace(0, np.pi / 3, 61):
    R = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]]); X = one @ R.T
    for s in np.arange(10, 40, 0.25):
        cands = [((i + 0.5 * j) * s, j * s * np.sqrt(3) / 2) for i in range(-6, 7) for j in range(-6, 7)]
        for off in ((0, 0), (s / 2, s * np.sqrt(3) / 4)):
            cs = [(cx + off[0], cy + off[1]) for cx, cy in cands]
            cs = [c for c in cs if np.sqrt(((X + c) ** 2).sum(1)).max() <= 49.5]
            if not cs: continue
            pts = [X + np.array(c) for c in cs]
            dmin = min(np.sqrt(((a[:, None] - b[None]) ** 2).sum(-1)).min() for a, b in it.combinations(pts, 2)) if len(pts) > 1 else 1e9
            if dmin < GAP: continue
            key = (len(cs), dmin)
            if best is None or key > best[0]: best = (key, rot, s, cs, X)
(nc, dmin), rot, s, cs, X = best
coords = {}
for k, c in enumerate(sorted(cs, key=lambda c: (round(c[1], 6), c[0]))):
    for i, f in enumerate(fam):
        coords[f'c{k}_{f}'] = (float(X[i, 0] + c[0]), float(X[i, 1] + c[1]))
allp = np.array(list(coords.values())); dd = np.sqrt(((allp[:, None] - allp[None]) ** 2).sum(-1)); np.fill_diagonal(dd, 1e9)
# PREFLIGHT on the tiled register (must be able to fail)
assert dd.min() >= 5.0, 'atoms closer than 5 um'
assert np.sqrt((allp ** 2).sum(1)).max() <= 50.0, 'outside 50 um radius'
assert len(coords) <= 100
inter = min(np.linalg.norm(np.array(coords[f'c{a}_{f}']) - np.array(coords[f'c{b}_{g}'])) for a, b in it.combinations(range(nc), 2) for f in fam for g in fam)
assert inter >= GAP
pre = dict(P['preflight']); pre.update({'copies': nc, 'atoms': len(coords), 'min_intercopy_um': inter, 'V_intercopy_max': C6 / inter ** 6, 'tiling': f'hex, rot {rot:.4f} rad, pitch {s} um'})
assert pre['V_intercopy_max'] < pre["deltaf"] / 40, 'inter-copy tail not negligible'
json.dump({'register_um': coords, 'preflight': pre, 'emulation': P['emulation']}, open(sys.argv[2], 'w'), indent=1)
print('copies', nc, 'atoms', len(coords), 'min intercopy', round(inter, 2), 'um; max radius', round(np.sqrt((allp ** 2).sum(1)).max(), 2), 'min any', round(dd.min(), 2))
