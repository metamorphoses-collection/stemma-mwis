"""EXPLORATORY (14 Nov 2026, not pre-registered): which single witness, if it were contaminated, explains the conflict
left by the exact rooted MWIS? Same families and thresholds as bedier_top.py (MINW 2, analysis witnesses of the frozen run).

Conflict mass M = total weight of families (weight >= CROSS_W = 5) that cross the MWIS. For each witness X, X is deleted
from every group and the MWIS is re-solved: drop_X = 1 - M_X / M. Deleting a witness that is cited often removes mass
anyway, so drop_X is set against involvement_X, the share of all family weight whose families contain X:
excess_X = drop_X / involvement_X (> 1: X carries more of the conflict than its share of the evidence).
Stability: B unit-bootstrap replicates; share of replicates in which each witness has the largest excess.
Usage: contam_witness.py text_id ... --out FILE [--B 100]"""
import json, os, argparse, random, collections, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rooted_mwis import mwis
B_ = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'historic', 'bedier')
MINW, CROSS_W = 2, 5
cross = lambda a, b: bool(a & b) and not (a <= b or b <= a)


def fams(U, W):
    F = []
    for u in U:
        fs = [frozenset(x for x in g if x in W) for g in u['groups'][1:]]
        fs = [g for g in fs if len(g) >= 2 and g != W]
        if fs: F.append(fs)
    return F


def mass(F, W):
    c = collections.Counter(f for fs in F for f in fs); Wt = {f: w for f, w in c.items() if w >= MINW}
    sel, _, st = mwis(Wt); assert st == 'OPTIMAL', st
    exc = {f: w for f, w in Wt.items() if w >= CROSS_W and f not in sel and any(cross(f, s) for s in sel)}
    return sum(exc.values()), sel, exc, Wt


def drop_witness(F, x):
    out = []
    for fs in F:
        g = [frozenset(f - {x}) for f in fs]; g = [f for f in g if len(f) >= 2]
        if g: out.append(g)
    return out


def analyse(F, W):
    M0, sel, exc, Wt = mass(F, W); tot = sum(Wt.values()); r = {}
    for x in sorted(W):
        Wx = W - {x}; Fx = [[f for f in fs if f != Wx] for fs in drop_witness(F, x)]; Fx = [fs for fs in Fx if fs]
        Mx = mass(Fx, Wx)[0]; inv = sum(w for f, w in Wt.items() if x in f) / tot
        d = 1 - Mx / M0 if M0 else 0.0
        r[x] = dict(M_without=Mx, drop=d, involvement=inv, excess=d / inv if inv else 0.0)
    return M0, sel, exc, r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('texts', nargs='+'); ap.add_argument('--out', required=True)
    ap.add_argument('--B', type=int, default=100); a = ap.parse_args(); res = {}
    key = lambda f, order: ''.join(sorted(f, key=order.index))
    for t in a.texts:
        d = json.load(open(os.path.join(B_, 'classification_a10', t + '.json'))); order = d['analysis_witnesses']; W = frozenset(order)
        F = fams(json.load(open(os.path.join(B_, 'runs_a10', t, 'units.json'))), W)
        M0, sel, exc, r = analyse(F, W)
        rng = random.Random(1); win = collections.Counter()
        for _ in range(a.B):
            Fb = [F[rng.randrange(len(F))] for _ in F]; _, _, _, rb = analyse(Fb, W)
            win[max(rb, key=lambda x: rb[x]['excess'])] += 1
        res[t] = dict(conflict_mass=M0, mwis=[key(f, order) for f in sorted(sel, key=lambda f: -len(f))],
                      excluded_crossing=sorted(([key(f, order), w] for f, w in exc.items()), key=lambda z: -z[1])[:12],
                      per_witness=r, top_suspect_bootstrap={x: v / a.B for x, v in win.most_common()})
        print(f'== {t}: conflict mass {M0}; strongest excluded', res[t]['excluded_crossing'][:6])
        for x in sorted(r, key=lambda x: -r[x]['excess']):
            print(f"   {x}: drop {r[x]['drop']:.2f}  involvement {r[x]['involvement']:.2f}  excess {r[x]['excess']:.2f}")
        print('   top suspect (bootstrap):', res[t]['top_suspect_bootstrap'])
    json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
