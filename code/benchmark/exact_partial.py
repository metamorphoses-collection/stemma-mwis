"""EXACT maximum-weight jointly completable partial splits (Theo, Theorem A.1), solved by cutting planes.

Master: MWIS on the overlap-conflict graph (a relaxation: every jointly completable family is overlap-compatible)
        + no-good cuts  sum_{i in K} y_i <= |K| - 1  for infeasible cores K.
Check:  can the chosen partial splits be completed together into pairwise compatible full splits?
        CP-SAT feasibility on the chosen splits only. Splits are oriented away from a reference witness r
        (the most observed one): cluster C_i = side without r; o_i = which side r falls on when r is unobserved.
        Compatible <=> C_i, C_j disjoint or nested (laminar). One assumption literal per split gives an
        infeasible core when the check fails.
Loop until the master's choice is completable: that choice is then optimal for the exact problem, and its
completions are the full splits (no nearest-neighbour extension, no greedy filter).
Weights: each split counted at most once per unit (no binary double count).

Usage: exact_partial.py <matrix.num> [--truth correct.txt --out dir] [--json out.json] [--tl 600]
"""
import sys, os, json, time, argparse, itertools
from collections import Counter
from ortools.sat.python import cp_model
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stemma_mwis import load_num, tree_from_clusters, write_graph


def partial_splits(rows, n):
    cnt = Counter()
    for r in rows:
        obs = [i for i in range(n) if r[i] != '?']; O = 0
        for i in obs: O |= 1 << i
        keys = set()
        for val in set(r[i] for i in obs):
            A = 0
            for i in obs:
                if r[i] == val: A |= 1 << i
            B = O & ~A
            if bin(A).count('1') >= 2 and bin(B).count('1') >= 2:
                keys.add((min(A, B), max(A, B)))
        for k in keys: cnt[k] += 1
    return cnt


def overlap_conflict(s, t):
    A1, B1 = s; A2, B2 = t
    return bool(A1 & A2) and bool(A1 & B2) and bool(B1 & A2) and bool(B1 & B2)


def master(nodes, W, edges, cuts, hint, tl):
    m = cp_model.CpModel(); y = [m.NewBoolVar('') for _ in nodes]
    for i, j in edges: m.AddBoolOr([y[i].Not(), y[j].Not()])
    for K in cuts: m.Add(sum(y[i] for i in K) <= len(K) - 1)
    m.Maximize(sum(w * v for w, v in zip(W, y)))
    if hint is not None:
        for i, v in enumerate(y): m.AddHint(v, i in hint)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl
    # deterministic by default (one worker, fixed seed): the reported optimum is chosen without reference to any truth
    s.parameters.num_workers = int(os.environ.get('MWIS_WORKERS', '1')); s.parameters.random_seed = int(os.environ.get('MWIS_SEED', '0'))
    st = s.Solve(m)
    return [i for i in range(len(nodes)) if s.Value(y[i])], s.StatusName(st), s.ObjectiveValue(), s.BestObjectiveBound()


def completion_check(nodes, S, n, r, tl, sim=None):
    """Feasibility of joint completion for splits S. Returns (True, clusters) or (False, core indices)."""
    m = cp_model.CpModel()
    a = {i: m.NewBoolVar(f'a{i}') for i in S}
    c = {}   # c[i][x]: Python bool (fixed) or BoolVar; True = x in the cluster of i (side without r)
    for i in S:
        A, B = nodes[i]; O = A | B
        o = bool(A >> r & 1) if O >> r & 1 else m.NewBoolVar('')   # True: r on side A, cluster = B side
        row = {}
        for x in range(n):
            if x == r: continue
            if A >> x & 1: row[x] = (not o) if isinstance(o, bool) else o.Not()
            elif B >> x & 1: row[x] = o
            else: row[x] = m.NewBoolVar('')
        c[i] = row
    def val(l):
        return int(l) if isinstance(l, bool) else None
    def neg(l):
        return (not l) if isinstance(l, bool) else l.Not()
    npairs = nclauses = 0
    for i, j in itertools.combinations(S, 2):
        Oi = nodes[i][0] | nodes[i][1]; Oj = nodes[j][0] | nodes[j][1]
        full = (1 << n) - 1
        if Oi == full and Oj == full: continue          # both complete: overlap test already decided it
        opts = []
        for kind in ('D', 'IJ', 'JI'):
            clauses = []; possible = True
            for x in range(n):
                if x == r: continue
                p, q = c[i][x], c[j][x]
                if kind == 'D': lits = [neg(p), neg(q)]
                elif kind == 'IJ': lits = [neg(p), q]
                else: lits = [p, neg(q)]
                if any(val(l) == 1 for l in lits): continue
                lits = [l for l in lits if val(l) is None]
                if not lits: possible = False; break
                clauses.append(lits)
            if not possible: continue
            if not clauses: opts = None; break               # this relation holds whatever the completion
            v = m.NewBoolVar('')
            for cl in clauses: m.AddBoolOr(cl).OnlyEnforceIf(v); nclauses += 1
            opts.append(v)
        if opts is None: continue
        npairs += 1
        if not opts: m.AddBoolOr([a[i].Not(), a[j].Not()])
        else: m.AddBoolOr(opts).OnlyEnforceIf([a[i], a[j]])
    if sim is not None:   # tie-break among valid completions: an unobserved witness sides with its most similar observed one
        terms = []
        for i in S:
            A, B = nodes[i]; O = A | B; obs = [x for x in range(n) if O >> x & 1]
            for x in range(n):
                if x == r or O >> x & 1: continue
                nb = max(obs, key=lambda j: sim[x][j])
                if nb == r: continue
                tgt = c[i][nb]
                eq = m.NewBoolVar('')
                if isinstance(tgt, bool):
                    m.Add(c[i][x] == int(tgt)).OnlyEnforceIf(eq)
                else:
                    m.Add(c[i][x] == tgt).OnlyEnforceIf(eq)
                terms.append(eq)
        m.Maximize(sum(terms))
        for i in S: m.Add(a[i] == 1)
    else:
        m.AddAssumptions([a[i] for i in S])
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 1; s.parameters.random_seed = 0
    st = s.Solve(m)
    info = dict(pairs=npairs, clauses=nclauses, status=s.StatusName(st))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        clusters = []
        for i in S:
            C = frozenset(x for x in range(n) if x != r and (val(c[i][x]) == 1 or (val(c[i][x]) is None and s.BooleanValue(c[i][x]))))
            clusters.append((i, C))
        return True, clusters, info
    if st == cp_model.INFEASIBLE:
        idx = s.SufficientAssumptionsForInfeasibility()
        lit2i = {a[i].Index(): i for i in S}
        core = [lit2i[k] for k in idx if k in lit2i] or list(S)
        return False, core, info
    raise RuntimeError(f'completion check undecided: {info}')


def solve(wit, rows, tl=600, log=print):
    global SIM
    from stemma_mwis import impute
    n = len(wit); _, SIM = impute(wit, rows)
    obsc = [sum(r[i] != '?' for r in rows) for i in range(n)]
    r = max(range(n), key=lambda i: obsc[i])
    cnt = partial_splits(rows, n); nodes = list(cnt); W = [cnt[k] for k in nodes]
    edges = [(i, j) for i, j in itertools.combinations(range(len(nodes)), 2) if overlap_conflict(nodes[i], nodes[j])]
    log(f'reference {wit[r]}  split types {len(nodes)}  instances {sum(W)}  overlap conflicts {len(edges)}')
    cuts, hint, it, t0 = [], None, 0, time.time()
    first_obj = None
    while True:
        it += 1
        S, st, obj, bd = master(nodes, W, edges, cuts, hint, tl)
        if first_obj is None: first_obj = obj
        ok, res, info = completion_check(nodes, S, n, r, tl)
        log(f'iter {it}: master {st} obj {obj:.0f} bound {bd:.0f} chosen {len(S)} | check {info} -> {"COMPLETABLE" if ok else f"core {len(res)}"}  ({time.time()-t0:.0f}s)')
        if ok:
            ok2, res2, info2 = completion_check(nodes, S, n, r, tl, sim=SIM)
            log(f'  similarity tie-break: {info2}')
            if ok2: res = res2
            return dict(reference=wit[r], split_types=len(nodes), instances=sum(W), overlap_conflicts=len(edges),
                        iterations=it, cuts=len(cuts), relaxation_obj=first_obj, exact_obj=obj, master_status=st,
                        chosen=len(S)), [C for _, C in res], r
        cuts.append(res); hint = set(S) - set(res[:1])


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('num'); ap.add_argument('--truth'); ap.add_argument('--out')
    ap.add_argument('--json'); ap.add_argument('--tl', type=int, default=600)
    a = ap.parse_args()
    wit, rows = load_num(a.num)
    summary, clusters, r = solve(wit, rows, a.tl)
    n = len(wit)
    # clusters are oriented away from r; tree_from_clusters expects clusters on the side without witness 0
    canon = []
    for C in clusters:
        if 0 in C: C = frozenset(range(n)) - C
        if 2 <= len(C) <= n - 2: canon.append(C)
    canon = list(dict.fromkeys(canon))
    summary['families'] = len(canon)
    if a.truth and a.out:
        from score import score
        os.makedirs(a.out, exist_ok=True)
        f = os.path.join(a.out, 'exact_partial.txt')
        write_graph(tree_from_clusters(canon, n), wit, f, directed=False)
        summary['skeleton'], summary['orientation'] = score(a.truth, f)
        print(f"skeleton {summary['skeleton']:.1f}%")
    summary['family_sets'] = [sorted(wit[i] for i in C) for C in canon]
    print(json.dumps({k: v for k, v in summary.items() if k != 'family_sets'}))
    if a.json: json.dump(summary, open(a.json, 'w'), indent=1)
