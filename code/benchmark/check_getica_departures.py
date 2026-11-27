"""Mommsen's Getica, joint departures of his orders and why the rooted optimum has no II+III family.

Definition (as in the paper's count, code/audit/getica_order_pairs.py in the MC repo, reimplemented here):
  a unit (edited reading identified, lemma_group == 0) is a JOINT DEPARTURE of orders P and Q when ONE departure
  group (a non-lemma reading) contains every PRESENT member of P and every present member of Q, and no present
  member of the third order; P and Q each have >= 1 present member. Order I = HPVL (A set aside), II = OB + S
  where present, III = XYZ. 'Present' = cited in some group of the unit (negative apparatus; Y, Z, S are known only
  where cited).
Then: the exact departure groups behind the II+III count, their rooted weights (rooted_mwis.py coding: a family is
the exact set of witnesses of a departure group, weight >= 2 units), and which optimum family each crosses.
Usage: check_getica_departures.py data/mommsen_getica/units.json data/mommsen_getica/rooted_lemma.json data/mommsen_getica/rooted_majority.json out.json
"""
import json, sys, collections, itertools
from ortools.sat.python import cp_model

U = json.load(open(sys.argv[1])); RL = json.load(open(sys.argv[2])); RM = json.load(open(sys.argv[3])); out = sys.argv[4]
known = [u for u in U if u['lemma_group'] == 0]
O = {'I': set('HPVL'), 'II': set('SOB'), 'III': set('XYZ')}
pairs = [('I', 'II', 'III'), ('I', 'III', 'II'), ('II', 'III', 'I')]
cnt = collections.Counter(); groups = {p[0] + '+' + p[1]: collections.Counter() for p in pairs}
for u in known:
    present = set(''.join(u['groups'])); P = {k: v & present for k, v in O.items()}
    hit = set()
    for g in u['groups'][1:]:
        g = set(g)
        for a, b, z in pairs:
            if P[a] and P[b] and P[a] | P[b] <= g and not (g & P[z]):
                k = a + '+' + b
                if k not in hit: hit.add(k); cnt[k] += 1      # count units, not groups
                groups[k][''.join(sorted(g, key='HPVLASOBXYZ'.index))] += 1
order = 'HPVLASOBXYZ'
key = lambda f: ''.join(sorted(f, key=order.index))
# rooted weights (lemma) as in rooted_mwis.families/weigh
F = []
for u in known:
    F.append([frozenset(g) for g in u['groups'][1:] if len(g) >= 2])
W = collections.Counter(f for fs in F for f in fs)
opt = [frozenset(r[0]) for r in RL['mwis']]
cross = lambda a, b: bool(a & b) and not (a <= b or b <= a)
II_III = []
for g, c in groups['II+III'].most_common():
    f = frozenset(g)
    II_III.append(dict(group=g, units=c, rooted_weight=W.get(f, 0),
                       crosses_optimum=[key(o) for o in opt if cross(f, o)]))
# what would it cost to force a single II+III family into the rooted optimum (lemma coding, weights >= 2)?
def rooted_opt(force=None):
    fam = [f for f, w in W.items() if w >= 2]
    if force is not None and force not in fam: fam.append(force)
    m = cp_model.CpModel(); x = {f: m.NewBoolVar('') for f in fam}
    for a, b in itertools.combinations(fam, 2):
        if cross(a, b): m.Add(x[a] + x[b] <= 1)
    if force is not None: m.Add(x[force] == 1)
    m.Maximize(sum(W.get(f, 0) * x[f] for f in fam))
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.random_seed = 0; s.parameters.max_time_in_seconds = 120
    st = s.Solve(m); return s.ObjectiveValue(), s.StatusName(st), [key(f) for f in fam if s.Value(x[f])]
base = rooted_opt()
forced = {}
for g in ['OBXYZ', 'SOBXYZ', 'OBXY', 'AOBXY', 'AOBXYZ', 'ASOBXYZ']:
    o, st, sel = rooted_opt(frozenset(g))
    forced[g] = dict(weight=o, status=st, loss=base[0] - o, kept=sel)
res = dict(units=len(U), edited_reading_identified=len(known), unidentified=len(U) - len(known),
           unidentified_examples=[dict(lemma=u['lemma'], groups=u['groups']) for u in U if u['lemma_group'] is None][:24],
           joint_departure_units=dict(cnt), departure_groups=dict((k, dict(v.most_common())) for k, v in groups.items()),
           II_III_groups=II_III, rooted_lemma_optimum=dict(weight=base[0], status=base[1], families=base[2]),
           forcing_a_II_III_family=forced,
           majority_reference=dict(units=RM['units'], candidate_families=RM['candidate_families'], weight=RM['weight'],
                                   status=RM['status'], families=RM['mwis'], stability_ge_0_1=RM['stability']),
           lemma_reference=dict(units=RL['units'], candidate_families=RL['candidate_families'], weight=RL['weight'],
                                families=RL['mwis']))
json.dump(res, open(out, 'w'), indent=1)
print('joint departure units:', dict(cnt))
print('II+III groups (group, units, rooted weight, crosses):')
for r in II_III[:25]: print('  ', r['group'], r['units'], r['rooted_weight'], r['crosses_optimum'])
print('I+II groups:', dict(groups['I+II'])); print('I+III groups:', dict(groups['I+III']))
print('rooted optimum (lemma):', base[0], base[1])
for g, v in forced.items(): print('force', g, '-> weight', v['weight'], 'loss', v['loss'])
