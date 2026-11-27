"""v2 of partial.py (11 Nov 2026): adds --dedupe (each split counted once per unit; fixes the binary double count found by
Theo) and --variants-only (drop lemma-vs-mixed-rest splits). partial.py itself is kept frozen for the 10 Nov results.
MWIS on PARTIAL splits: conflict tested only on witnesses observed in both loci (no imputation in selection)."""
import sys, itertools, json
from collections import Counter
from ortools.sat.python import cp_model
sys.path.insert(0, '.')
from stemma_mwis import load_num, impute, tree_from_clusters, write_graph
from score import score

def partial_splits(wit, rows, dedupe=False, varonly=False):
    n = len(wit); cnt = Counter()
    for r in rows:
        obs = [i for i in range(n) if r[i] != '?']
        O = sum(1 << i for i in obs)
        keys = []
        for val in set(r[i] for i in obs):
            A = sum(1 << i for i in obs if r[i] == val); B = O & ~A
            if bin(A).count('1') >= 2 and bin(B).count('1') >= 2:
                if varonly and len(set(r[i] for i in obs if r[i] != val)) > 1 and val == majority(r, obs): continue
                keys.append((min(A, B), max(A, B)))
        for k in (set(keys) if dedupe else keys): cnt[k] += 1
    return cnt

def majority(r, obs):
    from collections import Counter
    return Counter(r[i] for i in obs).most_common(1)[0][0]

def conflict(s, t):
    A1, B1 = s; A2, B2 = t
    return bool(A1 & A2) and bool(A1 & B2) and bool(B1 & A2) and bool(B1 & B2)

def solve(nodes, W, tl=300):
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in nodes]; E = 0
    for i, j in itertools.combinations(range(len(nodes)), 2):
        if conflict(nodes[i], nodes[j]):
            m.AddBoolOr([x[i].Not(), x[j].Not()]); E += 1
    m.Maximize(sum(w * v for w, v in zip(W, x)))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 8
    st = s.Solve(m)
    return [i for i in range(len(nodes)) if s.Value(x[i])], s.StatusName(st), s.ObjectiveValue(), s.BestObjectiveBound(), E

if __name__ == '__main__':
    num, truth, out = sys.argv[1:4]; dd = '--dedupe' in sys.argv; vo = '--variants-only' in sys.argv
    wit, raw = load_num(num); n = len(wit)
    _, sim = impute(wit, raw)
    cnt = partial_splits(wit, raw, dd, vo); nodes = list(cnt); W = [cnt[k] for k in nodes]
    ch, st, obj, bd, E = solve(nodes, W)
    print(f'partial splits {len(nodes)} instances {sum(W)} conflicts {E} status {st} obj {obj} bound {bd} chosen {len(ch)}')
    # extend chosen partial splits to full: each unobserved witness joins the side of its most similar observed witness
    full = []
    for i in sorted(ch, key=lambda k: -W[k]):
        A, B = nodes[i]; side = set(j for j in range(n) if A >> j & 1)
        obs = [j for j in range(n) if (A | B) >> j & 1]
        for m_ in range(n):
            if not (A | B) >> m_ & 1:
                nb = max(obs, key=lambda j: sim[m_][j])
                if A >> nb & 1: side.add(m_)
        side = frozenset(side)
        if 0 in side: side = frozenset(range(n)) - side
        if 2 <= len(side) <= n - 2: full.append((side, W[i]))
    kept = []
    for c, w in full:
        if all(not (c & d) or c <= d or d <= c for d, _ in kept): kept.append((c, w))
    print(f'extended {len(full)}, kept laminar {len(kept)} (dropped {len(full) - len(kept)}, weight {sum(w for _, w in full) - sum(w for _, w in kept)})')
    p = tree_from_clusters([c for c, _ in kept], n)
    f = f'{out}/partial_mwis.txt'; write_graph(p, wit, f, directed=False)
    print('skeleton %.1f%%' % score(truth, f)[0])
