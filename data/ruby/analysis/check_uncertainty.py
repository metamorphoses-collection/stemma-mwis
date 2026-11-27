"""Q3/Q7: uncertainty of the hardware shares, respecting that the copies of one shot are correlated.
For each run: pooled share of maximum sets and of conflict-free outcomes; naive binomial SE (copy samples treated as
independent); cluster-robust SE with the shot as cluster (equal cluster sizes: SE = sd of per-shot means / sqrt(shots));
design effect and intra-shot correlation; whole-shot bootstrap 95% interval; per-copy-position shares (min-max, chi-square
homogeneity); per pass for Priscus. Lawrence also: E, R and E/(E+R) with shot-clustered SE (delta method + shot bootstrap).
Output: checks_20261123/uncertainty.json"""
import math, json, itertools as it
import numpy as np
from collections import Counter, defaultdict
from check_common import GRAPHS, LAW_E, LAW_R, conflict_fn, maximum_sets, load_shots, dump

rng = np.random.default_rng(20261123); B = 20000


def chi2_sf(x, k):
    # survival of chi-square(k) via regularized gamma (series), adequate for reporting
    from math import lgamma, exp, log
    a, xx = k / 2, x / 2
    if xx <= 0: return 1.0
    s, term, n = 1 / a, 1 / a, 1
    while term > 1e-15 * s and n < 10000:
        term *= xx / (a + n); s += term; n += 1
    return max(0.0, 1 - exp(-xx + a * log(xx) - lgamma(a)) * s)


def summarise(per_shot, ncop):
    """per_shot: array (shots, ncop) of 0/1 indicators."""
    S, k = per_shot.shape
    m = per_shot.mean(1); p = per_shot.mean()
    naive = math.sqrt(p * (1 - p) / (S * k))
    clus = m.std(ddof=1) / math.sqrt(S)
    deff = (clus / naive) ** 2 if naive > 0 else float('nan')
    icc = (deff - 1) / (k - 1) if k > 1 else float('nan')
    ix = rng.integers(0, S, (B, S)); bm = m[ix].mean(1)
    cp = per_shot.mean(0)
    # chi-square homogeneity across copy positions (treats shots as independent within a column)
    obs = np.stack([per_shot.sum(0), S - per_shot.sum(0)]); exp_ = obs.sum(1, keepdims=True) * S / (S * k)
    chi = float(((obs - exp_) ** 2 / np.where(exp_ > 0, exp_, 1)).sum())
    return dict(share=p, se_naive_binomial=naive, se_shot_clustered=clus, design_effect=deff, intra_shot_corr=icc,
                ci95_shot_clustered=[p - 1.96 * clus, p + 1.96 * clus],
                ci95_shot_bootstrap=list(np.percentile(bm, [2.5, 97.5])),
                per_copy=[float(x) for x in cp], per_copy_min=float(cp.min()), per_copy_max=float(cp.max()),
                per_copy_sd=float(cp.std(ddof=1)), per_copy_sd_expected_if_homogeneous=math.sqrt(p * (1 - p) / S),
                chi2_copy_homogeneity=chi, chi2_df=k - 1, chi2_p=chi2_sf(chi, k - 1))


out = {}
for name, g in GRAPHS.items():
    conf = conflict_fn(g); mis = set(maximum_sets(g))
    copies, rows = load_shots(name)
    valid = lambda s: not any(conf(a, b) for a, b in it.combinations(s, 2))
    shots = []; files = []
    for f, outs, n in rows:
        for _ in range(n): shots.append(outs); files.append(f)
    k = len(copies); S = len(shots)
    M = np.array([[o in mis for o in sh] for sh in shots], dtype=float)
    V = np.array([[valid(o) for o in sh] for sh in shots], dtype=float)
    rec = dict(shots=S, copies=k, copy_samples=S * k, maximum=summarise(M, k), conflict_free=summarise(V, k))
    if name == 'priscus':
        fa = np.array(files); rec['per_pass'] = {}
        for f in g['files']:
            sel = fa == f
            rec['per_pass'][f] = dict(shots=int(sel.sum()), maximum=float(M[sel].mean()), conflict_free=float(V[sel].mean()),
                                      maximum_se_clustered=float(M[sel].mean(1).std(ddof=1) / math.sqrt(sel.sum())))
        deboor = frozenset({'BMP', 'MP'})
        D = np.array([[o == deboor for o in sh] for sh in shots], dtype=float)
        rec['deBoor_tree'] = summarise(D, k)
        rec['per_MIS'] = {}
        for s in sorted(mis, key=sorted):
            X = np.array([[o == s for o in sh] for sh in shots], dtype=float)
            r = summarise(X, k); rec['per_MIS']['{' + ','.join(x for x in g['fam'] if x in s) + '}'] = {kk: r[kk] for kk in ('share', 'se_shot_clustered', 'per_copy_min', 'per_copy_max')}
    if name == 'lai':
        ed = frozenset({'CG', 'DE', 'AB', 'DEF'})
        rec['per_MIS'] = {}
        for s in sorted(mis, key=sorted):
            X = np.array([[o == s for o in sh] for sh in shots], dtype=float)
            r = summarise(X, k); rec['per_MIS']['{' + ','.join(x for x in g['fam'] if x in s) + '}'] = {kk: r[kk] for kk in ('share', 'se_shot_clustered', 'per_copy_min', 'per_copy_max')}
        # tree vs runner-up, shot-clustered
        X1 = np.array([[o == ed for o in sh] for sh in shots], dtype=float)
        X2 = np.array([[o == frozenset({'CG', 'AB', 'DEF', 'DF'}) for o in sh] for sh in shots], dtype=float)
        d = (X1 - X2).mean(1)
        rec['tree_minus_runnerup'] = dict(diff=float(d.mean()), se_shot_clustered=float(d.std(ddof=1) / math.sqrt(S)))
    if name == 'lawrence':
        E = np.array([[o == LAW_E for o in sh] for sh in shots], dtype=float)
        R = np.array([[o == LAW_R for o in sh] for sh in shots], dtype=float)
        e, r = E.sum(1), R.sum(1)
        share = e.sum() / (e.sum() + r.sum())
        # delta-method ratio estimator, shot clusters
        mu_e, mu_r = e.mean(), r.mean(); t = mu_e + mu_r
        z = (e * mu_r - r * mu_e) / t ** 2
        se_ratio = z.std(ddof=1) / math.sqrt(S)
        ix = rng.integers(0, S, (B, S)); bs = e[ix].sum(1) / (e[ix].sum(1) + r[ix].sum(1))
        pc = [float(E[:, j].sum() / (E[:, j].sum() + R[:, j].sum())) for j in range(k)]
        rec['E'] = summarise(E, k); rec['R'] = summarise(R, k)
        rec['E_over_E_plus_R'] = dict(share=share, n_E=int(e.sum()), n_R=int(r.sum()), se_shot_clustered_delta=se_ratio,
                                      ci95_delta=[share - 1.96 * se_ratio, share + 1.96 * se_ratio],
                                      ci95_shot_bootstrap=list(np.percentile(bs, [2.5, 97.5])),
                                      per_copy=pc, per_copy_min=min(pc), per_copy_max=max(pc))
    out[name] = rec
    m, v = rec['maximum'], rec['conflict_free']
    print(f"{name}: shots {S} x {k} copies | maximum {m['share']:.4f} se naive {m['se_naive_binomial']:.4f} clustered {m['se_shot_clustered']:.4f}"
          f" (deff {m['design_effect']:.2f}, icc {m['intra_shot_corr']:.3f}) per-copy {m['per_copy_min']:.3f}-{m['per_copy_max']:.3f} chi2 p {m['chi2_p']:.2g}"
          f" | conflict-free {v['share']:.4f} se clustered {v['se_shot_clustered']:.4f} (deff {v['design_effect']:.2f}) per-copy {v['per_copy_min']:.3f}-{v['per_copy_max']:.3f} chi2 p {v['chi2_p']:.2g}")
    if 'per_pass' in rec: print('  passes:', {f[-8:-5]: (rec['per_pass'][f]['shots'], round(rec['per_pass'][f]['maximum'], 3), round(rec['per_pass'][f]['conflict_free'], 3)) for f in rec['per_pass']})
    if 'per_MIS' in rec: print('  per MIS:', {kk: (round(vv['share'], 4), round(vv['se_shot_clustered'], 4), round(vv['per_copy_min'], 3), round(vv['per_copy_max'], 3)) for kk, vv in rec['per_MIS'].items()})
    if 'tree_minus_runnerup' in rec: print('  Lai tree - runner-up:', rec['tree_minus_runnerup'])
    if 'E_over_E_plus_R' in rec:
        q = rec['E_over_E_plus_R']; print(f"  E {rec['E']['share']:.4f} (se cl {rec['E']['se_shot_clustered']:.4f}, copies {rec['E']['per_copy_min']:.3f}-{rec['E']['per_copy_max']:.3f}),"
                                         f" R {rec['R']['share']:.4f} (se cl {rec['R']['se_shot_clustered']:.4f}); E/(E+R) {q['share']:.4f} se {q['se_shot_clustered_delta']:.4f}"
                                         f" boot {q['ci95_shot_bootstrap'][0]:.4f}-{q['ci95_shot_bootstrap'][1]:.4f}; per copy {q['per_copy_min']:.3f}-{q['per_copy_max']:.3f}")
print('->', dump(out, 'uncertainty.json'))
