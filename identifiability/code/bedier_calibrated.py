r"""Coincidence-calibrated Bedier test on birth-death traditions (Proposition polytomynull applied).

A pruning rule at a nominal cutoff, evaluated empirically (no test level is claimed). For each family C of the raw
optimum with exactly two children D_1, D_2 in the estimated rooted tree, parent P and siblings, the second-order null
"C is a false resolution of a polytomy" (Proposition polytomynull) predicts the expected support
   mu_0 = theta q(D_1) q(D_2) L + (kappa - theta) q(P) q(S) L      (second term only if C has a single sibling S and
                                                                     P is not the root; theta = alpha^2 kappa)
with q(D) = w(D)/L the frequency with which D is exactly a non-archetype reading class (singletons included).
C is removed if P(Poisson(mu_0) >= w(C)) >= LEVEL. Families with more than two children are not tested (kept): no
null is proved for them. Removing families keeps the family laminar. theta is estimated from the collation: over leaf pairs (i,j) whose smallest
common family F in the raw tree has at least 4 witnesses and F minus {i,j} is not a family, theta_hat =
(sum C_ij + 1/2) / sum w_i w_j / L, with C_ij units where {i,j} is a reading class; kappa_hat likewise from R_ij,
units in which exactly i and j depart. Oracle arm: theta = alpha^2 kappa_true, kappa_true = 1 + 1/SHAPE (1 if equal).
Usage: bedier_calibrated.py OUT.json [N]; env as bedier_on_birthdeath.py."""
import sys, os, json, math
import numpy as np
from scipy.stats import poisson
src = open(os.path.join(os.path.dirname(__file__), 'bedier_on_birthdeath.py')).read().split('results = {}')[0]
OUTP = sys.argv[1]; NT = int(sys.argv[2]) if len(sys.argv) > 2 else 300
sys.argv = [sys.argv[0], '/dev/null', str(NT)]
exec(src)
LEVEL = 0.01


def unit_stats(T, surv, hot_filter=False):
    """w for every cluster of size 1..X-1, and per pair C_ij (class) and R_ij (exactly i, j depart).
    hot_filter: drop units with three or more distinct departing readings (third-order events in the small-rate
    regime; in practice the places that change many times). Returns also the number of units kept."""
    S = T[surv]; X = len(surv); w = {}; Cp = {}; Rp = {}; kept = 0
    for u in range(L):
        col = S[:, u]; dep = np.flatnonzero(col != 0)
        if hot_filter and len(set(col[dep].tolist())) >= 3: continue
        kept += 1
        for v in set(col[dep].tolist()):
            Cl = frozenset(np.flatnonzero(col == v).tolist())
            if 1 <= len(Cl) <= X - 1: w[Cl] = w.get(Cl, 0) + 1
        if len(dep) == 2:
            p = (int(dep[0]), int(dep[1])); Rp[p] = Rp.get(p, 0) + 1
            if col[dep[0]] == col[dep[1]]: Cp[p] = Cp.get(p, 0) + 1
    return w, Cp, Rp, kept


def tree_of(S, X):
    fam = sorted(S, key=len); full = frozenset(range(X))
    def parent(C): return min((F for F in fam if C < F), key=len, default=full)
    def children(C):
        subs = [F for F in fam if F < C]; maxs = [F for F in subs if not any(F < G for G in subs)]
        cov = set().union(*maxs) if maxs else set()
        return maxs + [frozenset([x]) for x in C - cov]
    return parent, children, full


def estimate(S, X, w, Cp, Rp):
    parent, children, full = tree_of(S, X); num = rn = den = 0.0
    for i in range(X):
        for j in range(i + 1, X):
            pair = frozenset((i, j)); F = min((G for G in S if pair <= G), key=len, default=full)
            if len(F) < 4 or (F - pair) in S: continue
            num += Cp.get((i, j), 0); rn += Rp.get((i, j), 0)
            den += w.get(frozenset([i]), 0) * w.get(frozenset([j]), 0) / LL
    if den == 0: return None, None
    return (num + 0.5) / den, (rn + 0.5) / den


def calibrated(S, X, w, theta, kappa):
    parent, children, full = tree_of(S, X); keep = set()
    for C in S:
        D = children(C); P = parent(C); sib = [F for F in children(P) if F != C]
        if len(D) != 2: keep.add(C); continue
        mu = theta * w.get(D[0], 0) * w.get(D[1], 0) / LL
        if len(sib) == 1 and P != full: mu += max(kappa - theta, 0.0) * w.get(P, 0) * w.get(sib[0], 0) / LL
        if poisson.sf(w.get(C, 0) - 1, max(mu, 1e-12)) < LEVEL: keep.add(C)
    return keep


kappa_true = 1 + 1 / SHAPE if SHAPE > 0 else 1.0
results = {}
for alpha in ALPHAS:
    arms = ['raw', 'calibrated', 'oracle', 'k>=3', 'filtered raw', 'filtered calibrated', 'filtered oracle']
    fallback = fallback_f = 0; th_f = []; kept_frac = []
    rows = {a: [] for a in ['true'] + arms}; acc = {a: [0, 0, 0] for a in arms}; th = []
    got = tries = 0
    while got < NT and tries < 200000:
        tries += 1; sim = simulate_tree()
        if sim is None: continue
        parent_, alive = sim; ns = sum(alive)
        if not 3 <= ns <= 30: continue
        surv, truth = prune(parent_, alive); T = texts(parent_, alpha)
        w1, Cp, Rp, LL = unit_stats(T, surv)
        wc = {C: x for C, x in w1.items() if 2 <= len(C) <= ns - 1}
        raw = mwis(wc)
        t_hat, k_hat = estimate(raw, ns, w1, Cp, Rp)
        if t_hat is None: t_hat, k_hat = 1.0, 1.0; fallback += 1
        th.append(t_hat)
        rec = {'raw': raw, 'calibrated': calibrated(raw, ns, w1, t_hat, k_hat),
               'oracle': calibrated(raw, ns, w1, alpha ** 2 * kappa_true, kappa_true),
               'k>=3': mwis({C: x for C, x in wc.items() if x >= 3})}
        globals()['LL'] = LL
        w2, Cp2, Rp2, LL = unit_stats(T, surv, hot_filter=True); globals()['LL'] = LL
        wc2 = {C: x for C, x in w2.items() if 2 <= len(C) <= ns - 1}
        raw2 = mwis(wc2); t2, k2 = estimate(raw2, ns, w2, Cp2, Rp2)
        if t2 is None: t2, k2 = 1.0, 1.0; fallback_f += 1
        th_f.append(t2); kept_frac.append(LL / L)
        rec['filtered raw'] = raw2; rec['filtered calibrated'] = calibrated(raw2, ns, w2, t2, k2)
        rec['filtered oracle'] = calibrated(raw2, ns, w2, alpha ** 2 * kappa_true, kappa_true)
        rows['true'].append(metrics(truth, ns))
        for a, S in rec.items():
            rows[a].append(metrics(S, ns)); acc[a][0] += len(S & truth); acc[a][1] += len(S - truth); acc[a][2] += len(truth)
        got += 1
    results[f'alpha={alpha}'] = dict({a: summarise(v) for a, v in rows.items()},
                                     recovery={a: dict(true_found=x[0], false_added=x[1], true_total=x[2]) for a, x in acc.items()},
                                     theta_true=alpha ** 2 * kappa_true, theta_hat_median=float(np.median(th)),
                                     theta_hat_iqr=[float(np.percentile(th, 25)), float(np.percentile(th, 75))],
                                     theta_hat_filtered_median=float(np.median(th_f)),
                                     theta_hat_filtered_iqr=[float(np.percentile(th_f, 25)), float(np.percentile(th_f, 75))],
                                     fallback_unfiltered=fallback, fallback_filtered=fallback_f,
                                     units_kept_after_filter_median=float(np.median(kept_frac)))
    r = results[f'alpha={alpha}']
    print(f'alpha={alpha} theta {r["theta_true"]:.2f} hat {r["theta_hat_median"]:.2f}', {a: (r[a]['bifid_root'], r[a]['i3'], round(r['recovery'][a]['true_found'] / r['recovery'][a]['true_total'], 3), round(r['recovery'][a]['false_added'] / NT, 2)) for a in ['true'] + arms if a in r['recovery'] or a == 'true'} if False else '', flush=True)
    for a in ['true'] + arms:
        x = r[a]; rc = r['recovery'].get(a)
        print(f"  {a:10} bifid {x['bifid_root']:.2f} deg2 {x['p_deg2']:.2f} deg3 {x['p_deg3']:.2f} i3 {x['i3']:.2f}",
              f"| true {rc['true_found'] / rc['true_total']:.3f} false/trad {rc['false_added'] / NT:.2f}" if rc else '', flush=True)
results['params'] = dict(L=L, e=E, gamma_shape=SHAPE, level=LEVEL, traditions_per_alpha=NT, kappa_true=kappa_true)
json.dump(results, open(OUTP, 'w'), indent=1)
