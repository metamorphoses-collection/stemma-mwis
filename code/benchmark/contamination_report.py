"""Contamination report: what the exact MWIS rejected, and which witnesses explain the conflict.

For every candidate family (split) NOT kept by the exact solution, with weight >= --min-weight:
  - the kept families it conflicts with (on the completed tree);
  - the smallest set of witnesses (size 1, then 2) whose removal makes it compatible with every kept family.
Those witnesses are the candidates for contamination (or for a misplaced family): a philologist's work list.
Also per witness: how much rejected weight its removal would rescue.
Usage: contamination_report.py <matrix.num> [--min-weight 3] [--json out.json]
"""
import sys, os, json, itertools, argparse
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from stemma_mwis import load_num
import exact_partial

ap = argparse.ArgumentParser(); ap.add_argument('num'); ap.add_argument('--min-weight', type=int, default=3)
ap.add_argument('--json'); ap.add_argument('--tl', type=int, default=300)
a = ap.parse_args()
wit, rows = load_num(a.num); n = len(wit)
_, kept_clusters, r = exact_partial.solve(wit, rows, a.tl, log=lambda *x: None)
FULL = frozenset(range(n))
kept = []
for C in kept_clusters:
    C = frozenset(C)
    if r in C: C = FULL - C
    if 2 <= len(C) <= n - 2 and C not in kept: kept.append(C)

cnt = exact_partial.partial_splits(rows, n)


def conflicts(A, B, C, drop=frozenset()):
    """partial split A|B (observed) vs kept full cluster C (side without r), ignoring witnesses in drop."""
    A2, B2 = A - drop, B - drop
    if len(A2) < 1 or len(B2) < 1: return False
    C2 = C - drop; D2 = (FULL - C) - drop
    return bool(A2 & C2) and bool(A2 & D2) and bool(B2 & C2) and bool(B2 & D2)


bits = lambda m: frozenset(i for i in range(n) if m >> i & 1)
kept_set = set(kept)
report, rescue = [], Counter()
for (Am, Bm), w in sorted(cnt.items(), key=lambda x: -x[1]):
    if w < a.min_weight: break
    A, B = bits(Am), bits(Bm)
    clash = [C for C in kept if conflicts(A, B, C)]
    if not clash: continue                       # compatible with the tree: kept or implied
    expl = None
    for k in (1, 2):
        for D in itertools.combinations(sorted(A | B), k):
            D = frozenset(D)
            if all(not conflicts(A, B, C, D) for C in clash):
                expl = [sorted(wit[i] for i in D)] if expl is None else expl + [sorted(wit[i] for i in D)]
        if expl: break
    side = A if r not in A else B
    report.append(dict(family=sorted(wit[i] for i in side), against=sorted(wit[i] for i in (A | B) - side),
                       weight=w, conflicts_with=[sorted(wit[i] for i in C) for C in clash][:6],
                       n_conflicts=len(clash), explained_by_removing=expl or []))
    for e in (expl or []):
        if len(e) == 1: rescue[e[0]] += w
print(f'kept families {len(kept)}; rejected families with weight >= {a.min_weight}: {len(report)}')
for x in report[:25]:
    print(f"w={x['weight']:3d}  {'+'.join(x['family'])}  vs {x['n_conflicts']} kept  -> remove {x['explained_by_removing'][:4] or 'none (k<=2)'}")
print('\nwitnesses whose single removal rescues most rejected weight:', rescue.most_common(12))
if a.json: json.dump(dict(kept=[sorted(wit[i] for i in C) for C in kept], rejected=report, rescue=rescue.most_common()), open(a.json, 'w'), indent=1)
