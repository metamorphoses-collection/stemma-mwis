"""Review A17 + A25 (Donatus): margins of the editor's 12 families under the exact unrooted method, the editor's
stemma weight against the optimum, and tree-likeness indices.

m(sigma) = W_forb(sigma) - W_imp(sigma), where
  W_forb: exact optimum of Eq. (1) with every candidate vertex that PROPOSES sigma removed (a vertex proposes sigma
          when, on the witnesses it observes, one of its sides equals sigma; these are the units counted as
          "support" in Table 2);
  W_imp:  exact optimum with sigma imposed as a full split (a fully observed vertex of weight 0 forced in; its
          supporters stay available, every candidate crossing it is excluded, joint completion includes it).
For a family of the optimum, the margin by which it wins is W_imp - W_forb (= -m). Sign convention as asked.
Also, for every family of the reported optimum, the same test (forbid its supporters): if W_forb < W_opt the
family is supported in EVERY optimum.
Usage: check_donatus_margins.py <donatus.num> <out.json> [--workers 8]
Outputs sigla and counts only (no readings).
"""
import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'donatus'))
from check_exact_lib import Instance
from build_matrix import GROUPS

ap = argparse.ArgumentParser(); ap.add_argument('num'); ap.add_argument('out')
ap.add_argument('--workers', type=int, default=8); ap.add_argument('--tl', type=int, default=600)
a = ap.parse_args()
t0 = time.time()
I = Instance(a.num); n = I.n
print(f'{n} witnesses, {len(I.rows)} units, {len(I.nodes)} split types, {sum(I.W)} instances, {len(I.edges)} conflicts, ref {I.wit[I.r]}')
base = I.solve(workers=a.workers, log=print)
Wopt = base['obj']; fams = I.canon_families(base['clusters'])
print(f'optimum {Wopt} ({base["status"]}, bound {base["bound"]}), {len(fams)} families, {base["seconds"]:.1f}s')
names = lambda C: sorted(I.wit[i] for i in C)
famset = set(fams)

ed = {}
for g, mem in GROUPS.items():
    Cset = frozenset(I.wit.index(w) for w in mem)
    C = I.mask(mem)
    sup = I.supporters(C); crs = I.crossers(C)
    f = I.solve(forbid=sup, workers=a.workers)
    im = I.solve(impose=[C], workers=a.workers, tiebreak=False)
    fam_f = set(I.canon_families(f['clusters']))
    ed[g] = dict(size=len(mem), support=sum(I.W[i] for i in sup), support_vertices=len(sup),
                 crossing=sum(I.W[i] for i in crs), in_optimum=Cset in famset,
                 W_forb=f['obj'], W_imp=im['obj'], margin_forb_minus_imp=f['obj'] - im['obj'],
                 status=(f['status'], im['status']), still_in_tree_when_forbidden=Cset in fam_f)
    print(f"{g:8s} n={len(mem):2d} sup={ed[g]['support']:3d} cross={ed[g]['crossing']:3d} in_opt={ed[g]['in_optimum']!s:5} "
          f"W_forb={f['obj']:.0f} W_imp={im['obj']:.0f} m={f['obj']-im['obj']:+.0f}", flush=True)

# editor's stemma as a whole
masks = {g: I.mask(m) for g, m in GROUPS.items()}
sets = {g: frozenset(I.wit.index(w) for w in m) for g, m in GROUPS.items()}
compat = all(not (sets[x] & sets[y]) or sets[x] <= sets[y] or sets[y] <= sets[x] for x in sets for y in sets)
allimp = I.solve(impose=list(masks.values()), workers=a.workers, tiebreak=False)
# weight of the evidence the editor's tree accommodates: candidates compatible with all 12 families
crossing_any = set()
for C in masks.values(): crossing_any |= set(I.crossers(C))
compatible_weight = sum(w for i, w in enumerate(I.W) if i not in crossing_any)
subsets = {'all_but_nu': [g for g in GROUPS if g != 'nu'],
           'all_but_nu_upsilon': [g for g in GROUPS if g not in ('nu', 'upsilon')],
           'pi_and_tau': ['pi', 'tau']}
sub = {}
for k, gs in subsets.items():
    r_ = I.solve(impose=[masks[g] for g in gs], workers=a.workers, tiebreak=False)
    sub[k] = dict(families=gs, W_imp=r_['obj'], status=r_['status'], loss_vs_opt=Wopt - r_['obj'])
    print('impose', k, r_['obj'], flush=True)
editor = dict(pairwise_compatible=compat, sum_of_supports=sum(v['support'] for v in ed.values()),
              W_opt=Wopt, W_imp_all12=allimp['obj'], status_all12=allimp['status'], loss_all12=Wopt - allimp['obj'],
              W_imp_subsets=sub, weight_of_candidates_compatible_with_all12=compatible_weight)
print('editor:', editor)

# every family of the reported optimum: is it supported in every optimum?
persist = []
for C in fams:
    M = sum(1 << i for i in C)
    sup = I.supporters(M)
    f = I.solve(forbid=sup, workers=a.workers, tiebreak=False)
    persist.append(dict(family=names(C), support=sum(I.W[i] for i in sup), W_forb=f['obj'], loss=Wopt - f['obj'],
                        in_every_optimum=f['obj'] < Wopt))
print('families of the optimum supported in every optimum:', sum(p['in_every_optimum'] for p in persist), 'of', len(persist))

# A25 tree-likeness
total = sum(I.W)
inf, cons = I.consistent_units(fams); inf2, cons2 = I.consistent_units_vertexwise(base['S'])
tl = dict(units=len(I.rows), split_types=len(I.nodes), total_candidate_weight=total, W_opt=Wopt,
          share_of_candidate_weight=Wopt / total, informative_units=inf, units_consistent_with_tree=cons,
          share_consistent=cons / inf, units_consistent_vertexwise=cons2)
print('tree-likeness:', tl)
json.dump(dict(base=dict(obj=Wopt, status=base['status'], bound=base['bound'], relaxation=base['relaxation'],
                         iterations=base['iterations'], cuts=base['cuts'], seconds=base['seconds'], families=len(fams),
                         family_sets=[names(C) for C in fams]),
               editor_families=ed, editor_stemma=editor, optimum_family_persistence=persist, treelikeness=tl,
               workers=a.workers, seconds=time.time() - t0), open(a.out, 'w'), indent=1)
