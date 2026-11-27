"""Calibrated encoding of family weights as per-atom final detunings (used by calib_calibrate.py and calib_run.py).

Goal: Gibbs log-odds between maximal histories h, k equal theta * (W_h - W_k), theta = beta_eff * s, i.e. final detuning
difference s per unit of support. Rules (rad/us):
 * families present in EVERY maximal history ('universal') cancel from all shares: detuning delta_f (max yield);
 * non-universal family i: delta_i = s * w_i + lambda * o_i, where o >= 1 is a SIZE-NEUTRAL offset vector
   (sum_{i in h} o_i equal for every maximal history h, found by LP) so the lambda term adds the same reward to every
   history and cancels from the shares. lambda puts the smallest non-universal detuning at the floor delta_min.
   (A plain per-family floor a would add a*delta_f*(|h|-|k|) to the log-odds: a size bias. It is harmless only when all
   competing histories have the same size: Priscus, Lai. Lawrence (4 vs 3 families) needs o = (B+R 1, P+R 2, B+F+R 1).)
 * if the largest non-universal detuning would exceed delta_f, s is reduced until it fits (saturation; recorded).
"""
import numpy as np
from scipy.optimize import linprog
from calib_common import INST, DELTAF, maximal_sets

DMIN = 0.3 * DELTAF   # floor = the preps' a = 0.3 level, 1.70 rad/us


def neutral_offsets(name):
    g = INST[name]; F = g['fam']; mx = maximal_sets(g)
    uni = [f for f in F if all(f in s for s in mx)]; J = [f for f in F if f not in uni]
    M = np.array([[f in s for f in J] for s in mx], float)
    # variables o_J (>=1), c ; M o - c = 0 ; minimise sum o
    nJ = len(J); A_eq = np.hstack([M, -np.ones((len(mx), 1))]); b_eq = np.zeros(len(mx))
    r = linprog(np.r_[np.ones(nJ), 0], A_eq=A_eq, b_eq=b_eq, bounds=[(1, None)] * nJ + [(0, None)], method='highs')
    if not r.success: return uni, J, np.ones(nJ), False
    return uni, J, r.x[:nJ], True


def encode(name, w, s, dmin=DMIN, dmax=DELTAF):
    """w: dict family -> support; s: requested rad/us per unit support. Returns (w_hat vector in fam order, s_used, info)."""
    g = INST[name]; F = g['fam']; uni, J, o, ok = neutral_offsets(name)
    wj = np.array([w[f] for f in J], float)
    def build(s):
        lam = max((dmin - s * wj) / o); return s * wj + lam * o
    s_used = s; d = build(s)
    if d.max() > dmax + 1e-9:
        lo, hi = 0.0, s
        for _ in range(60):
            mid = (lo + hi) / 2
            if build(mid).max() <= dmax: lo = mid
            else: hi = mid
        s_used = lo; d = build(lo)
    delta = {f: dmax for f in uni}; delta.update({f: float(x) for f, x in zip(J, d)})
    what = np.array([delta[f] / DELTAF for f in F])
    return what, s_used, dict(universal=uni, offsets=dict(zip(J, o.tolist())), offsets_size_neutral=ok,
                              delta={f: delta[f] for f in F}, saturated=bool(s_used < s - 1e-12))


def raw_encoding(name, w, a=0.3):
    """Control: the preps' encoding w_eff = a + (1-a) w / w_max on every family (here with the case's own weights)."""
    F = INST[name]['fam']; m = max(w.values())
    return np.array([a + (1 - a) * w[f] / m for f in F])


def history_tails(name, P):
    from calib_common import interactions
    g = INST[name]; F = g['fam']; V = interactions(P); mx = maximal_sets(g)
    return np.array([sum(V[i, j] for i in range(len(F)) for j in range(i + 1, len(F)) if F[i] in s and F[j] in s) for s in mx])


def encode_tc(name, w, s, P, dmin=DMIN, dmax=DELTAF):
    """TAIL-COMPENSATED encoding (used for all calibrated runs). The compatible pairs of the checked layouts sit at
    >= 1.25 R_b, where V = C6/r^6 is still up to ~2.5 rad/us, so the true end-of-sweep reward of a history,
    eps_h = sum_{i in h} delta_i - T_h (T_h = vdW tails inside h), differs from the nominal one by up to ~3 rad/us
    between histories -- as large as the calibrated differences themselves. Solve instead, by LP,
        sum_{i in h} delta_i - T_h = s * W_h + c     for every maximal history h,   dmin <= delta_i <= dmax,
    maximising sum delta (yield). The common c cancels from the shares (this also absorbs the size-neutral floor).
    If infeasible, s is reduced to the largest feasible value (saturation)."""
    g = INST[name]; F = g['fam']; mx = maximal_sets(g); n = len(F)
    M = np.array([[f in h for f in F] for h in mx], float); W = M @ np.array([w[f] for f in F], float); Th = history_tails(name, P)
    def solve(s):
        A_eq = np.hstack([M, -np.ones((len(mx), 1))]); b_eq = s * W + Th
        r = linprog(np.r_[-np.ones(n), 0], A_eq=A_eq, b_eq=b_eq, bounds=[(dmin, dmax)] * n + [(None, None)], method='highs')
        return r.x[:n] if r.success else None
    s_used = s; d = solve(s)
    if d is None:
        lo, hi = 0.0, s
        if solve(0.0) is None: raise RuntimeError(f'{name}: tails cannot be compensated within [{dmin},{dmax}]')
        for _ in range(50):
            mid = (lo + hi) / 2
            if solve(mid) is not None: lo = mid
            else: hi = mid
        s_used = lo; d = solve(lo)
    return d / DELTAF, s_used, dict(delta=dict(zip(F, d.tolist())), saturated=bool(s_used < s - 1e-9), tails=Th.tolist())


def encode_qp(name, w, s, P, bias=None, beta=None, dmin=DMIN, dmax=DELTAF):
    """As encode_tc, but (i) the detuning vector is the one CLOSEST to all-at-delta_f (least squares, so it moves
    continuously with s instead of jumping between LP vertices -- an LP vertex jump at s = 0 flipped the Lawrence
    calibration point by ~2.7 in logit, i.e. the device responds to individual detunings, not only to the history
    rewards), and (ii) optional ZERO-POINT correction: bias[h] = log device share of history h measured on the same register
    with s = 0 (all history rewards equal); the targets become eps_h = s W_h - bias_h / beta + c.
    Saturation: s reduced to the largest feasible value (LP feasibility), recorded."""
    from scipy.optimize import minimize
    g = INST[name]; F = g['fam']; mx = maximal_sets(g); n = len(F)
    M = np.array([[f in h for f in F] for h in mx], float); W = M @ np.array([w[f] for f in F], float); Th = history_tails(name, P)
    b = np.zeros(len(mx)) if bias is None else np.asarray(bias, float) / beta
    A_eq = np.hstack([M, -np.ones((len(mx), 1))])
    def feas(s):
        r = linprog(np.zeros(n + 1), A_eq=A_eq, b_eq=s * W - b + Th, bounds=[(dmin, dmax)] * n + [(None, None)], method='highs')
        return r.success
    s_used = s; frac = 1.0
    if not feas(0.0):   # zero point not fully reachable within the detuning window: apply the largest feasible fraction
        b0 = b.copy(); lo, hi = 0.0, 1.0
        for _ in range(50):
            mid = (lo + hi) / 2; b = mid * b0
            if feas(0.0): lo = mid
            else: hi = mid
        frac = lo; b = lo * b0
    if not feas(s):
        lo, hi = 0.0, s
        for _ in range(50):
            mid = (lo + hi) / 2
            if feas(mid): lo = mid
            else: hi = mid
        s_used = lo
    rhs = s_used * W - b + Th
    r0 = linprog(np.r_[-np.ones(n), 0], A_eq=A_eq, b_eq=rhs, bounds=[(dmin, dmax)] * n + [(None, None)], method='highs')
    x0 = r0.x
    q = minimize(lambda x: float(((x[:n] - dmax) ** 2).sum()), x0, jac=lambda x: np.r_[2 * (x[:n] - dmax), 0],
                 constraints=[dict(type='eq', fun=lambda x: A_eq @ x - rhs, jac=lambda x: A_eq)],
                 bounds=[(dmin, dmax)] * n + [(None, None)], method='SLSQP', options=dict(ftol=1e-12, maxiter=500))
    x = q.x if q.success and np.abs(A_eq @ q.x - rhs).max() < 1e-6 else x0
    d = np.clip(x[:n], dmin, dmax)
    return d / DELTAF, s_used, dict(delta=dict(zip(F, d.tolist())), saturated=bool(s_used < s - 1e-9), tails=Th.tolist(),
                                    qp_ok=bool(q.success), bias_rad=(b).tolist(), bias_fraction=frac)
