"""Rooted (Lachmannian) MWIS on an edited apparatus: families are groups of >= 2 witnesses sharing a departure
from the editor's text; two families conflict when they overlap without nesting. Weight = number of units.
Exact MWIS by CP-SAT, then B bootstrap replicates over units for stability.

Only units whose editor's-text group is known (lemma_group == 0) are used; shared correct readings are
ignored, as Lachmann's rule requires.
Usage: rooted_mwis.py units.json --out DIR [--B 200] [--min 2] [--drop S,..]
"""
import json, argparse, random, os, collections
from ortools.sat.python import cp_model


def families(U, drop=()):
    F = []
    for u in U:
        if u.get('lemma_group') != 0: continue
        fs = []
        for g in u['groups'][1:]:
            g = ''.join(w for w in g if w not in drop)
            if len(g) >= 2: fs.append(frozenset(g))
        F.append(fs)
    return F


def mwis(weights):
    fam = [f for f, w in weights.items() if w > 0]
    m = cp_model.CpModel(); x = {f: m.NewBoolVar(''.join(sorted(f))) for f in fam}
    for i, a in enumerate(fam):
        for b in fam[i + 1:]:
            if (a & b) and not (a <= b or b <= a): m.Add(x[a] + x[b] <= 1)
    m.Maximize(sum(weights[f] * x[f] for f in fam))
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.random_seed = 0
    s.parameters.max_time_in_seconds = 120
    st = s.Solve(m)
    return {f for f in fam if s.Value(x[f])}, s.ObjectiveValue(), s.StatusName(st)


def weigh(F, mn):
    c = collections.Counter(f for fs in F for f in fs)
    return {f: w for f, w in c.items() if w >= mn}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('units'); ap.add_argument('--out', required=True)
    ap.add_argument('--B', type=int, default=200); ap.add_argument('--min', type=int, default=2)
    ap.add_argument('--drop', default=''); ap.add_argument('--order', default='HPVLASOBXYZ')
    a = ap.parse_args(); order = a.order
    U = json.load(open(a.units)); F = families(U, set(filter(None, a.drop.split(','))))
    W = weigh(F, a.min)
    sel, obj, status = mwis(W)
    key = lambda f: ''.join(sorted(f, key=order.index))
    rng = random.Random(1); stab = collections.Counter()
    for _ in range(a.B):
        Fb = [F[rng.randrange(len(F))] for _ in F]
        sb, _, _ = mwis(weigh(Fb, a.min))
        stab.update(sb)
    res = dict(units=len(F), candidate_families=len(W), status=status, weight=obj,
               mwis=[[key(f), W[f], stab[f] / a.B] for f in sorted(sel, key=lambda f: -W[f])],
               stability=[[key(f), W.get(f, 0), v / a.B] for f, v in stab.most_common() if v / a.B >= 0.1],
               top_weights=[[key(f), w] for f, w in sorted(W.items(), key=lambda t: -t[1])[:30]])
    os.makedirs(a.out, exist_ok=True)
    json.dump(res, open(os.path.join(a.out, 'rooted.json'), 'w'), indent=1)
    print(f"{res['units']} units, {res['candidate_families']} candidate families, {status}, weight {obj}")
    print('MWIS (family, shared departures, stability):')
    for r in res['mwis']: print('  ', r)
    print('top weights:', res['top_weights'][:20])


if __name__ == '__main__':
    main()
