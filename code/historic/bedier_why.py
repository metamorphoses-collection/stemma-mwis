"""EXPLORATORY (decided after the Bédier classification of 2026-11-13; not pre-registered): why are tops two-way?

Check 1, merge: the data's top blocks are the maximal stable families of the frozen run plus witnesses in none (k_D).
Every shared departure (group of >= 2 analysis witnesses, not all of them) that touches >= 2 blocks is a 'join'. A join
inside one editor branch is evidence for that branch as a unit above its blocks; a join that meets both branches is
evidence against the editor's split. Under a true two-way top, inside-joins should clearly outnumber across-joins.

Check 2, choice: in units where one editor branch departs as a whole (a majority of its witnesses in departure groups)
while the other branch is intact (none of its witnesses cited), the editor's text follows the intact branch. Counts per
branch show whether the text leans on one branch (Bédier: the two-way stemma leaves every disagreement to the editor).

Usage: bedier_why.py text_id ... --out FILE   (reads classification_a10/<t>.json and runs_a10/<t>/units.json)
"""
import json, os, argparse, collections
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'historic', 'bedier')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('texts', nargs='+'); ap.add_argument('--out', required=True); a = ap.parse_args()
    res = {}
    for t in a.texts:
        d = json.load(open(os.path.join(B, 'classification_a10', t + '.json')))
        W = set(d['analysis_witnesses']); top = [frozenset(x) for x in d['editor_top']]
        stable = [frozenset(k) for k, w, s in d['stable_families']]
        maximal = [f for f in stable if not any(f < g for g in stable)]
        blocks = maximal + [frozenset([w]) for w in W - set().union(*maximal)] if maximal else [frozenset([w]) for w in W]
        U = json.load(open(os.path.join(B, 'runs_a10', t, 'units.json')))
        joins = collections.Counter(); choice = collections.Counter(); split_units = 0
        for u in U:
            groups = [frozenset(x for x in g if x in W) for g in u['groups'][1:]]
            groups = [g for g in groups if len(g) >= 2 and g != W]
            cited = set().union(*[frozenset(x for x in g if x in W) for g in u['groups'][1:]]) if u['groups'][1:] else set()
            for g in groups:
                touched = [b for b in blocks if b & g]
                if len(touched) < 2: continue
                inside = [i for i, br in enumerate(top) if g <= br]
                meets = [i for i, br in enumerate(top) if g & br]
                if inside: joins[f'inside_{"".join(sorted(top[inside[0]]))}'] += 1; joins['inside'] += 1
                elif len(meets) >= 2: joins['across'] += 1
                else: joins['outside_top'] += 1
            if len(top) == 2:
                dep = [len(cited & br) * 2 > len(br) for br in top]          # branch departs as a whole (majority)
                intact = [not (cited & br) for br in top]
                for i in (0, 1):
                    if dep[i] and intact[1 - i]:
                        choice[''.join(sorted(top[1 - i]))] += 1; split_units += 1
        # size-matched null for the merge check: inside share of joins under every split with the editor's branch sizes
        jg = []
        for u in U:
            for g in u['groups'][1:]:
                g = frozenset(x for x in g if x in W)
                if len(g) >= 2 and g != W and len([b for b in blocks if b & g]) >= 2: jg.append(g)
        null_p = None; null_mean = None; share = None
        if len(top) == 2 and joins.get('inside', 0) + joins.get('across', 0) > 0:
            import itertools
            share = joins.get('inside', 0) / (joins.get('inside', 0) + joins.get('across', 0))
            s1, s2 = sorted(len(x) for x in top); Wl = sorted(W); shares = []; seen = set()
            jc = collections.Counter(jg)
            from math import comb
            import random
            def score(A_, B_):
                i_ = x_ = 0
                for g, n in jc.items():
                    if g <= A_ or g <= B_: i_ += n
                    elif g & A_ and g & B_: x_ += n
                if i_ + x_: shares.append(i_ / (i_ + x_))
            if comb(len(Wl), s1) * comb(len(Wl) - s1, s2) <= 50000:          # exact enumeration
                for c1 in itertools.combinations(Wl, s1):
                    rest = [w for w in Wl if w not in c1]
                    for c2 in itertools.combinations(rest, s2):
                        key = frozenset([frozenset(c1), frozenset(c2)])
                        if key in seen: continue
                        seen.add(key); score(set(c1), set(c2))
            else:                                                               # Monte Carlo, 20,000 splits, seed 1
                rng = random.Random(1)
                for _ in range(20000):
                    p_ = rng.sample(Wl, s1 + s2); score(set(p_[:s1]), set(p_[s1:]))
            null_mean = sum(shares) / len(shares); null_p = sum(1 for v in shares if v >= share - 1e-12) / len(shares)
        multi = len(top) == 2 and min(len(x) for x in top) >= 2
        res[t] = dict(merge_inside_share=share, merge_null_mean=null_mean, merge_p=null_p, choice_informative=multi,
                      editor_top=d['editor_top'], classification=d['classification'], k_D=d['k_D'],
                      blocks=[''.join(sorted(b)) for b in blocks], joins=dict(joins),
                      text_follows=dict(choice), split_units=split_units)
        ins, acr = joins.get('inside', 0), joins.get('across', 0)
        ch = ' / '.join(f'{k} {v}' for k, v in choice.items())
        ms = '' if share is None else f"inside share {share:.2f} (random {null_mean:.2f}, p {null_p:.3f})"
        print(f"{t:30s} {d['classification']:13s} kD {d['k_D']:2d} joins in {ins:4d} across {acr:4d} {ms} | follows: {ch if multi else '(single-witness branch)'}")
    json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
