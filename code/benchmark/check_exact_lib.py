"""Review helpers (23 Nov 2026): the exact partial-data MWIS of exact_partial.py, with extra constraints.

Re-uses exact_partial.partial_splits / overlap_conflict / completion_check unchanged; the master is rebuilt here
so that vertices can be forbidden (y_v = 0), a full split can be imposed (added as a fully observed vertex of
weight 0 with y = 1), and no-good cuts on whole solutions can be added for enumeration.
Nothing in the existing scripts or outputs is modified.
"""
import os, sys, time, itertools
from ortools.sat.python import cp_model
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_partial import partial_splits, overlap_conflict, completion_check
from stemma_mwis import load_num, impute


class Instance:
    def __init__(self, path):
        self.wit, self.rows = load_num(path)
        self.n = n = len(self.wit)
        _, self.sim = impute(self.wit, self.rows)
        obsc = [sum(r[i] != '?' for r in self.rows) for i in range(n)]
        self.r = max(range(n), key=lambda i: obsc[i])            # same reference as exact_partial.solve
        cnt = partial_splits(self.rows, n)
        self.nodes = list(cnt); self.W = [cnt[k] for k in self.nodes]
        self.full = (1 << n) - 1
        self.edges = [(i, j) for i, j in itertools.combinations(range(len(self.nodes)), 2)
                      if overlap_conflict(self.nodes[i], self.nodes[j])]
        # per-unit keys (the partial splits each unit proposes), for consistency counts
        self.unit_keys = []
        for row in self.rows:
            obs = [i for i in range(n) if row[i] != '?']; O = sum(1 << i for i in obs); keys = set()
            for val in set(row[i] for i in obs):
                A = sum(1 << i for i in obs if row[i] == val); B = O & ~A
                if bin(A).count('1') >= 2 and bin(B).count('1') >= 2: keys.add((min(A, B), max(A, B)))
            self.unit_keys.append(keys)

    def mask(self, names):
        return sum(1 << self.wit.index(w) for w in names)

    def supporters(self, C):
        """Vertices that propose family C (C | rest): on their observed witnesses, one side equals C."""
        out = []
        for i, (A, B) in enumerate(self.nodes):
            O = A | B; c = C & O
            if c == A or c == B: out.append(i)
        return out

    def crossers(self, C):
        s = (C, self.full & ~C)
        return [i for i, v in enumerate(self.nodes) if overlap_conflict(v, s)]

    def solve(self, forbid=(), impose=(), extra_cuts=(), tl=600, workers=8, seed=0, log=None, max_iter=200, tiebreak=True):
        """Exact optimum (cutting planes) with vertices `forbid` excluded and full splits `impose` (masks) forced.
        Returns dict(obj, status, bound, S (chosen vertex ids, imposed excluded), clusters, iterations, cuts, seconds)."""
        nodes = list(self.nodes); W = list(self.W); edges = list(self.edges)
        imp_ids = []
        for C in impose:
            s = (min(C, self.full & ~C), max(C, self.full & ~C))
            k = len(nodes); nodes.append(s); W.append(0); imp_ids.append(k)
            edges += [(i, k) for i in range(len(self.nodes)) if overlap_conflict(self.nodes[i], s)]
            edges += [(j, k) for j in imp_ids[:-1] if overlap_conflict(nodes[j], s)]
        cuts = [list(K) for K in extra_cuts]; t0 = time.time(); it = 0; first = None
        while True:
            it += 1
            m = cp_model.CpModel(); y = [m.NewBoolVar('') for _ in nodes]
            for i, j in edges: m.AddBoolOr([y[i].Not(), y[j].Not()])
            for K in cuts: m.Add(sum(y[i] for i in K) <= len(K) - 1)
            for i in forbid: m.Add(y[i] == 0)
            for k in imp_ids: m.Add(y[k] == 1)
            m.Maximize(sum(w * v for w, v in zip(W, y)))
            s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl
            s.parameters.num_workers = workers; s.parameters.random_seed = seed
            st = s.Solve(m); stn = s.StatusName(st)
            if st != cp_model.OPTIMAL:   # only a proven master optimum certifies; never read values otherwise
                return dict(obj=None, status=stn, bound=s.BestObjectiveBound() if st == cp_model.FEASIBLE else None,
                            S=[], clusters=[], iterations=it, cuts=len(cuts), seconds=time.time() - t0, certified=False)
            S = [i for i in range(len(nodes)) if s.Value(y[i])]
            obj, bd = s.ObjectiveValue(), s.BestObjectiveBound()
            if first is None: first = obj
            ok, res, info = completion_check(nodes, S, self.n, self.r, tl)
            if log: log(f'  iter {it}: {stn} obj {obj:.0f} bound {bd:.0f} -> {"ok" if ok else "core %d" % len(res)} ({time.time()-t0:.0f}s)')
            if ok:
                if tiebreak:   # same similarity tie-break among valid completions as exact_partial.solve
                    ok2, res2, _ = completion_check(nodes, S, self.n, self.r, tl, sim=self.sim)
                    if ok2: res = res2
                return dict(obj=obj, status=stn, bound=bd, certified=True, relaxation=first, S=[i for i in S if i < len(self.nodes)],
                            clusters=[C for _, C in res], iterations=it, cuts=len(cuts) - len(extra_cuts),
                            new_cuts=cuts[len(extra_cuts):], seconds=time.time() - t0)
            cuts.append(res)
            if it >= max_iter: raise RuntimeError('too many iterations')

    def canon_families(self, clusters):
        """Clusters (oriented away from r) -> canonical families oriented away from witness 0, non-trivial."""
        n = self.n; out = []
        for C in clusters:
            C = frozenset(C)
            if 0 in C: C = frozenset(range(n)) - C
            if 2 <= len(C) <= n - 2: out.append(C)
        return list(dict.fromkeys(out))

    def consistent_units(self, families):
        """Informative units (>= 1 proposed partial split) none of whose proposed splits crosses a tree family."""
        splits = [(sum(1 << i for i in C), self.full & ~sum(1 << i for i in C)) for C in families]
        inf = cons = 0
        for keys in self.unit_keys:
            if not keys: continue
            inf += 1
            if all(not overlap_conflict(k, s) for k in keys for s in splits): cons += 1
        return inf, cons

    def consistent_units_vertexwise(self, S):
        """Completion-free variant: every proposed split of the unit is chosen or crosses no chosen vertex."""
        Sset = set(S); Sn = [self.nodes[i] for i in S]
        idx = {v: i for i, v in enumerate(self.nodes)}
        inf = cons = 0
        for keys in self.unit_keys:
            if not keys: continue
            inf += 1
            if all(idx[k] in Sset or not any(overlap_conflict(k, t) for t in Sn) for k in keys): cons += 1
        return inf, cons
