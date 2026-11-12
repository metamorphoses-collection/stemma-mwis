"""Prototype: stemma as exact maximum-weight independent set over a split-conflict graph.
Pipeline: matrix -> (impute missing by nearest witness) -> splits per locus -> merge identical
splits into weighted 'family types' -> conflict edges (four-gamete) -> exact MWIS (CP-SAT) ->
laminar clusters -> tree -> optional ancestor contraction -> CASC graph file.
Also a neighbour-joining baseline on the same imputed matrix, for calibration."""
import sys, itertools, statistics, json
from collections import Counter, defaultdict
from ortools.sat.python import cp_model


def load_num(path):
    L = [l for l in open(path).read().split('\n') if l.strip()]
    wit = L[0].split()
    rows = [l.split() for l in L[1:]]
    return wit, rows


def impute(wit, rows):
    """Fill '?' with the reading of the most similar witness that has one (agreement on co-observed loci)."""
    n = len(wit)
    sim = [[0.0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            co = [(r[a], r[b]) for r in rows if r[a] != '?' and r[b] != '?']
            sim[a][b] = sum(x == y for x, y in co) / len(co) if co else 0.0
    order = [sorted((b for b in range(n) if b != a), key=lambda b: -sim[a][b]) for a in range(n)]
    out = []
    for r in rows:
        rr = list(r)
        for a in range(n):
            if rr[a] == '?':
                for b in order[a]:
                    if r[b] != '?':
                        rr[a] = r[b]
                        break
        out.append(rr)
    return out, sim


def splits(wit, rows, ref=0):
    """One split per reading shared by >=2 witnesses with >=2 outside. Canonical side = side without witness[ref]."""
    n = len(wit)
    cnt = Counter()
    for r in rows:
        for val in set(r):
            side = frozenset(i for i in range(n) if r[i] == val)
            if 2 <= len(side) <= n - 2:
                if ref in side:
                    side = frozenset(range(n)) - side
                cnt[side] += 1
    return cnt


def compatible(a, b, n):
    # canonical sides (both exclude ref): compatible iff nested or disjoint
    return not (a & b) or a <= b or b <= a


def mwis(nodes, weights, n, time_limit=120):
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f'x{i}') for i in range(len(nodes))]
    E = 0
    for i, j in itertools.combinations(range(len(nodes)), 2):
        if not compatible(nodes[i], nodes[j], n):
            m.AddBoolOr([x[i].Not(), x[j].Not()])
            E += 1
    m.Maximize(sum(int(w) * xi for w, xi in zip(weights, x)))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = 8
    st = s.Solve(m)
    chosen = [i for i in range(len(nodes)) if s.Value(x[i])]
    return chosen, s.StatusName(st), s.ObjectiveValue(), s.BestObjectiveBound(), E


def tree_from_clusters(clusters, n):
    """Laminar clusters (+ full set) -> parent map. Hidden node per cluster; leaves = witnesses."""
    cl = sorted(set(clusters), key=len)
    full = frozenset(range(n))
    allc = cl + [full]
    parent = {}
    for i, c in enumerate(allc[:-1]):
        parent[('h', i)] = next(('h', j) if j < len(allc) - 1 else ('h', 'root')
                                for j in range(i + 1, len(allc)) if c < allc[j])
    for w in range(n):
        smallest = next((i for i, c in enumerate(allc[:-1]) if w in c), None)
        parent[('w', w)] = ('h', smallest) if smallest is not None else ('h', 'root')
    return parent


def private_rate(wit, raw):
    n = len(wit)
    out = {}
    for a in range(n):
        obs = [r for r in raw if r[a] != '?']
        priv = sum(1 for r in obs if sum(1 for x in r if x == r[a]) == 1)
        out[a] = priv / max(1, len(obs))
    return out


def contract(parent, wit, raw, factor=0.25):
    """Pre-declared rule: a witness whose private-reading rate is <= factor x median rate, and which is the
    lowest-rate leaf under its parent, replaces that hidden parent (= a surviving exemplar)."""
    pr = private_rate(wit, raw)
    med = statistics.median(pr.values())
    kids = defaultdict(list)
    for c, p in parent.items():
        kids[p].append(c)
    merged = {}
    for p, ch in kids.items():
        leaves = [c[1] for c in ch if c[0] == 'w']
        if not leaves:
            continue
        best = min(leaves, key=lambda w: pr[w])
        if pr[best] <= factor * med:
            merged[p] = ('w', best)
    newp = {}
    for c, p in parent.items():
        if c in merged.values() and parent[c] in merged and merged[parent[c]] == c:
            pp = parent.get(p)
            if pp is not None:
                newp[c] = merged.get(pp, pp)
            continue
        if c in merged:
            continue
        newp[c] = merged.get(p, p)
    for h, w in merged.items():
        if h in parent and w not in newp:
            pp = parent[h]
            newp[w] = merged.get(pp, pp)
    return newp, merged, pr, med


def write_graph(parent, wit, path, directed=True):
    nodes = set(parent) | set(parent.values())
    def lab(x):
        return wit[x[1]] if x[0] == 'w' else f'{900 + (x[1] if x[1] != "root" else 99)}'
    names = sorted({lab(x) for x in nodes})
    idx = {nm: i for i, nm in enumerate(names)}
    M = [[-1] * len(names) for _ in names]
    for i in range(len(names)):
        M[i][i] = 0
    for c, p in parent.items():
        a, b = idx[lab(p)], idx[lab(c)]
        M[a][b] = 1
        if not directed:
            M[b][a] = 1
    with open(path, 'w') as f:
        f.write(f'{len(names)}\n' + '\n'.join(names) + '\n')
        for row in M:
            f.write(' '.join(str(v) for v in row) + '\n')


def nj(wit, rows):
    n = len(wit)
    D = {(a, b): sum(r[a] != r[b] for r in rows) / len(rows) for a in range(n) for b in range(n)}
    active = [('w', i) for i in range(n)]
    dist = {(x, y): D[x[1], y[1]] for x in active for y in active}
    parent, k = {}, 0
    while len(active) > 2:
        m = len(active)
        R = {x: sum(dist[x, y] for y in active if y != x) for x in active}
        best = min(((x, y) for x, y in itertools.combinations(active, 2)),
                   key=lambda p: (m - 2) * dist[p] - R[p[0]] - R[p[1]])
        u = ('h', k); k += 1
        x, y = best
        parent[x] = u; parent[y] = u
        for z in active:
            if z not in best:
                d = 0.5 * (dist[x, z] + dist[y, z] - dist[x, y])
                dist[u, z] = dist[z, u] = d
        dist[u, u] = 0
        active = [z for z in active if z not in best] + [u]
    a, b = active
    parent[a] = b
    return parent


if __name__ == '__main__':
    num, out = sys.argv[1], sys.argv[2]
    wit, raw = load_num(num)
    rows, _ = impute(wit, raw)
    n = len(wit)
    cnt = splits(wit, rows)
    nodes = list(cnt)
    w = [cnt[s] for s in nodes]
    chosen, status, obj, bound, E = mwis(nodes, w, n)
    clusters = [nodes[i] for i in chosen]
    parent = tree_from_clusters(clusters, n)
    write_graph(parent, wit, f'{out}/mwis_tree.txt', directed=False)
    p2, merged, pr, med = contract(parent, wit, raw)
    write_graph(p2, wit, f'{out}/mwis_contracted.txt', directed=False)
    write_graph(nj(wit, rows), wit, f'{out}/nj.txt', directed=False)
    json.dump({'witnesses': n, 'loci': len(raw), 'distinct_splits': len(nodes),
               'split_instances': sum(w), 'conflict_edges': E, 'status': status,
               'objective': obj, 'bound': bound, 'chosen_splits': len(chosen),
               'chosen_weight_share': obj / sum(w),
               'contracted': {wit[v[1]]: str(h) for h, v in merged.items()},
               'median_private_rate': med},
              open(f'{out}/run.json', 'w'), indent=1)
    print(open(f'{out}/run.json').read())
