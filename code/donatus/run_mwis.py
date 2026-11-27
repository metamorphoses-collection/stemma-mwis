"""Exact MWIS (partial-data version) on the Donatus matrix, then comparison with Licciardello's stemma.
Usage: run_mwis.py <matrix.num> [--ref Gb] [--min-weight 1] [--json out.json]
Prints aggregates only (family sets and weights are sigla, no readings)."""
import sys, os, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'benchmark'))
from stemma_mwis import load_num, impute
from partial import partial_splits, solve
from build_matrix import GROUPS

ap = argparse.ArgumentParser()
ap.add_argument('num'); ap.add_argument('--ref', default='Gb')
ap.add_argument('--min-weight', type=int, default=1); ap.add_argument('--json')
ap.add_argument('--tl', type=int, default=600)
ap.add_argument('--dedupe', action='store_true', help='count each split at most once per unit (no binary double count)')
ap.add_argument('--variants-only', action='store_true', help='splits only from variant readings (r != 0) or from a homogeneous complement')
a = ap.parse_args()

wit, raw = load_num(a.num)
r = wit.index(a.ref)
order = [r] + [i for i in range(len(wit)) if i != r]  # reference first: clusters = side without ref
wit = [wit[i] for i in order]; raw = [[row[i] for i in order] for row in raw]
n = len(wit)
_, sim = impute(wit, raw)
from collections import Counter
def splits_v2(wit, rows, dedupe, varonly):
    n = len(wit); cnt = Counter()
    for r in rows:
        obs = [i for i in range(n) if r[i] != '?']; O = sum(1 << i for i in obs); keys = []
        for val in set(r[i] for i in obs):
            A = sum(1 << i for i in obs if r[i] == val); B = O & ~A
            if bin(A).count('1') < 2 or bin(B).count('1') < 2: continue
            if varonly and val == '0' and len(set(r[i] for i in obs if r[i] != val)) > 1: continue
            keys.append((min(A, B), max(A, B)))
        for k in (set(keys) if dedupe else keys): cnt[k] += 1
    return cnt
cnt = splits_v2(wit, raw, a.dedupe, a.variants_only) if (a.dedupe or a.variants_only) else partial_splits(wit, raw)
nodes = [k for k in cnt if cnt[k] >= a.min_weight]; W = [cnt[k] for k in nodes]
ch, st, obj, bd, E = solve(nodes, W, tl=a.tl)
print(f'loci {len(raw)}  split types {len(nodes)}  instances {sum(W)}  conflict edges {E}  status {st}  obj {obj:.0f}  bound {bd:.0f}  chosen {len(ch)}')

full = []
for i in sorted(ch, key=lambda k: -W[k]):
    A, B = nodes[i]; side = set(j for j in range(n) if A >> j & 1)
    obs = [j for j in range(n) if (A | B) >> j & 1]
    for m in range(n):
        if not (A | B) >> m & 1:
            nb = max(obs, key=lambda j: sim[m][j])
            if A >> nb & 1: side.add(m)
    side = frozenset(side)
    if 0 in side: side = frozenset(range(n)) - side
    if 2 <= len(side) <= n - 2: full.append((side, W[i]))
merged = {}
for c, w in full:  # partial splits extending to the same full family are one family
    merged[c] = merged.get(c, 0) + w
full = sorted(merged.items(), key=lambda x: -x[1])
kept = []
for c, w in full:
    if all(not (c & d) or c <= d or d <= c for d, _ in kept): kept.append((c, w))
dropped = len(full) - len(kept)
print(f'extended {len(full)}  kept laminar {len(kept)}  dropped {dropped}')

names = lambda s: sorted(wit[i] for i in s)
L = {g: frozenset(wit.index(w) for w in ms) for g, ms in GROUPS.items()}
keptsets = {c: w for c, w in kept}

print('\nLicciardello group -> recovered?  (best kept family by Jaccard)')
recall = 0
res_groups = {}
for g, s in sorted(L.items(), key=lambda x: -len(x[1])):
    best = max(kept, key=lambda cw: len(cw[0] & s) / len(cw[0] | s)) if kept else (frozenset(), 0)
    j = len(best[0] & s) / len(best[0] | s) if kept else 0
    exact = s in keptsets
    recall += exact
    extra = names(best[0] - s); missing = names(s - best[0])
    res_groups[g] = dict(size=len(s), exact=exact, jaccard=round(j, 3), weight=best[1], extra=extra, missing=missing)
    print(f'  {g:8s} n={len(s):2d}  {"EXACT" if exact else "     "}  J={j:.2f}  w={best[1]:3d}  +{extra} -{missing}')
print(f'recall (exact): {recall}/{len(L)}')

print('\nKept families (weight desc): his group / nested in / crosses')
prec = 0; fams = []
for c, w in sorted(kept, key=lambda x: -x[1]):
    tag = next((g for g, s in L.items() if s == c), None)
    if tag: kind = f'= {tag}'; prec += 1
    else:
        inside = [g for g, s in L.items() if c < s]
        if inside:
            kind = 'inside ' + min(inside, key=lambda g: len(L[g])); prec += 1
        else:
            unions = [g for g, s in L.items() if s < c]
            kind = 'CROSSES' + (f' (contains {",".join(sorted(unions))})' if unions else '')
    fams.append(dict(members=names(c), weight=w, relation=kind))
    print(f'  w={w:3d}  {kind:30s} {names(c)}')
print(f'precision (= or nested in one of his groups): {prec}/{len(kept)}')
if a.json:
    json.dump(dict(loci=len(raw), split_types=len(nodes), instances=sum(W), conflict_edges=E,
                   status=st, obj=obj, bound=bd, chosen=len(ch), kept=len(kept), dropped=dropped,
                   groups=res_groups, families=fams, recall=recall, precision=prec), open(a.json, 'w'), indent=1)
