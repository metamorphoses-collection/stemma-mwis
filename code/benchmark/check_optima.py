"""Review A8 + A25: solver status, bound, time, number of optimal solutions, families common to all optima, and
tree-likeness, for one tradition, with the exact method of exact_partial.py (cutting planes; CP-SAT).

1. One exact solve with the given workers/seed (status, objective, bound, iterations, cuts, wall time).
2. Enumeration of optimal vertex sets: after each optimum S, add the no-good cut sum_{i in S} y_i <= |S|-1 and
   re-solve (cores found on the way are kept); stop when the master's optimum falls below W* (enumeration
   complete) or at --cap optima / --budget seconds (enumeration incomplete: a lower bound on the count).
   Distinct vertex sets can give the same tree; distinct trees (canonical families after the similarity
   tie-break completion) are also counted, and skeleton scores are given when --truth is supplied.
3. Families common to all optima: for each family C of the first optimum, the supporting vertices of C are
   forbidden; if the exact optimum then falls below W*, every optimum contains a unit proposing C (C is
   'forced'); otherwise an equally heavy solution without C exists.
4. Tree-likeness: W* against the total candidate weight; informative units none of whose proposed splits
   crosses a family of the reported tree.
Usage: check_optima.py <matrix.num> <out.json> [--truth correct.txt] [--workers 8] [--cap 200] [--budget 3600]
       [--no-enum] [--no-forced]
"""
import sys, os, json, time, argparse, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_exact_lib import Instance
from stemma_mwis import tree_from_clusters, write_graph

ap = argparse.ArgumentParser(); ap.add_argument('num'); ap.add_argument('out'); ap.add_argument('--truth')
ap.add_argument('--workers', type=int, default=8); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--cap', type=int, default=200); ap.add_argument('--budget', type=float, default=3600)
ap.add_argument('--tl', type=int, default=600)
ap.add_argument('--no-enum', action='store_true'); ap.add_argument('--no-forced', action='store_true')
a = ap.parse_args()
log = lambda *x: print(*x, flush=True)
I = Instance(a.num); n = I.n; names = lambda C: sorted(I.wit[i] for i in C)
log(f'{n} witnesses, {len(I.rows)} units, {len(I.nodes)} split types, {sum(I.W)} instances, {len(I.edges)} conflicts, ref {I.wit[I.r]}')


def skel(fams):
    if not a.truth: return None
    from score import score
    f = tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False).name
    write_graph(tree_from_clusters(list(fams), n), I.wit, f, directed=False)
    s = score(a.truth, f)[0]; os.unlink(f); return round(s, 2)


base = I.solve(workers=a.workers, seed=a.seed, tl=a.tl, log=log)
Wst = base['obj']; fams = I.canon_families(base['clusters'])
res = dict(matrix=os.path.basename(a.num), witnesses=n, units=len(I.rows), split_types=len(I.nodes),
           total_candidate_weight=sum(I.W), conflicts=len(I.edges), reference=I.wit[I.r],
           solver=dict(engine='OR-Tools CP-SAT', workers=a.workers, seed=a.seed, time_limit_per_solve_s=a.tl),
           solve=dict(obj=Wst, relaxation=base['relaxation'], status=base['status'], bound=base['bound'],
                      gap=base['bound'] - Wst, iterations=base['iterations'], cuts=base['cuts'],
                      seconds=round(base['seconds'], 1), families=len(fams), skeleton=skel(fams),
                      family_sets=[names(C) for C in fams]))
log('solve:', {k: v for k, v in res['solve'].items() if k != 'family_sets'})
inf, cons = I.consistent_units(fams); inf2, cons2 = I.consistent_units_vertexwise(base['S'])
res['treelikeness'] = dict(W_opt=Wst, total_candidate_weight=sum(I.W), share=Wst / sum(I.W),
                           informative_units=inf, consistent_with_tree=cons, share_consistent=cons / inf,
                           consistent_vertexwise=cons2)
log('treelikeness:', res['treelikeness'])

if not a.no_enum:
    t0 = time.time(); sols = [base]; cuts = list(base.get('new_cuts', [])); complete = False; stop = ''
    while True:
        if len(sols) >= a.cap: stop = 'cap'; break
        if time.time() - t0 > a.budget: stop = 'budget'; break
        nog = [sorted(s['S']) for s in sols]
        r = I.solve(extra_cuts=cuts + nog, workers=a.workers, seed=a.seed, tl=a.tl)
        cuts += r.get('new_cuts', [])
        if r['obj'] is None or r['obj'] < Wst - 1e-9:
            complete = (r['status'] == 'OPTIMAL') if r['obj'] is not None else (r['status'] == 'INFEASIBLE')
            stop = f"next best {r['obj']} ({r['status']}, bound {r['bound']})"; break
        if r['status'] != 'OPTIMAL': stop = f"unproven at {r['obj']} ({r['status']})"; break
        sols.append(r)
        if len(sols) % 10 == 0: log(f'  optima found {len(sols)} ({time.time()-t0:.0f}s)')
    trees = {}
    for s in sols:
        key = frozenset(I.canon_families(s['clusters'])); trees.setdefault(key, 0); trees[key] += 1
    tscores = [skel(k) for k in trees] if a.truth else None
    vint = set.intersection(*[set(s['S']) for s in sols])
    tint = frozenset.intersection(*trees.keys())
    res['enumeration'] = dict(optimal_vertex_sets=len(sols), complete=complete, stopped=stop,
                              seconds=round(time.time() - t0, 1), distinct_trees=len(trees),
                              tree_scores=tscores, vertices_in_all_found=len(vint),
                              vertices_per_optimum=sorted(set(len(s['S']) for s in sols)),
                              families_in_all_found_trees=[names(C) for C in tint],
                              families_per_tree=sorted(set(len(k) for k in trees)))
    log('enumeration:', {k: v for k, v in res['enumeration'].items() if k != 'families_in_all_found_trees'})

if not a.no_forced:
    t0 = time.time(); out = []
    for C in fams:
        M = sum(1 << i for i in C); sup = I.supporters(M)
        f = I.solve(forbid=sup, workers=a.workers, seed=a.seed, tl=a.tl, tiebreak=False)
        out.append(dict(family=names(C), support=sum(I.W[i] for i in sup), W_forb=f['obj'], status=f['status'],
                        bound=f['bound'], forced=(f['bound'] is not None and f['bound'] < Wst - 1e-9) or
                        (f['obj'] is not None and f['obj'] < Wst - 1e-9 and f['status'] == 'OPTIMAL')))
        log(f"  {''.join(out[-1]['family'])[:40]:40s} sup {out[-1]['support']:3d} W_forb {f['obj']} {f['status']} -> {'FORCED' if out[-1]['forced'] else 'optional'}")
    res['families_forced'] = dict(forced=sum(o['forced'] for o in out), of=len(out), seconds=round(time.time() - t0, 1),
                                  detail=out)
json.dump(res, open(a.out, 'w'), indent=1)
