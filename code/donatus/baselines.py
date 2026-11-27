"""Baselines on the identical Donatus matrix, scored exactly like run_mwis.py against Licciardello's groups.
- nj:      neighbour-joining on p-distance over co-observed units (no imputation)
- newick:  any unrooted tree (parsimony from pars.R, etc.)
- edges:   an edge list (Camps & Cafiero PCC output: two columns of node names)
Every tree is converted to its clusters = clades on the side away from the reference witness (Gb).
Usage:
  baselines.py nj <matrix.num>
  baselines.py newick <matrix.num> <tree.nwk>
  baselines.py edges <matrix.num> <edges.tsv>
  add --json out.json to save the score."""
import sys, os, re, json, itertools
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmark'))
from stemma_mwis import load_num
from build_matrix import GROUPS

REF = 'Gb'


def nj_tree(wit, rows):
    n = len(wit)
    D = [[0.0] * n for _ in range(n)]
    for a, b in itertools.combinations(range(n), 2):
        co = [(r[a], r[b]) for r in rows if r[a] != '?' and r[b] != '?']
        D[a][b] = D[b][a] = sum(x != y for x, y in co) / len(co) if co else 0.5
    adj = defaultdict(set)
    active = list(range(n)); dist = {(i, j): D[i][j] for i in range(n) for j in range(n)}
    nxt = n
    while len(active) > 2:
        m = len(active)
        R = {x: sum(dist[x, y] for y in active if y != x) for x in active}
        x, y = min(itertools.combinations(active, 2), key=lambda p: (m - 2) * dist[p] - R[p[0]] - R[p[1]])
        u = nxt; nxt += 1
        adj[u] |= {x, y}; adj[x].add(u); adj[y].add(u)
        for z in active:
            if z not in (x, y):
                dist[u, z] = dist[z, u] = 0.5 * (dist[x, z] + dist[y, z] - dist[x, y])
        dist[u, u] = 0
        active = [z for z in active if z not in (x, y)] + [u]
    a, b = active; adj[a].add(b); adj[b].add(a)
    return adj


def parse_newick(s, wit):
    s = re.sub(r':[0-9.eE+-]+', '', s.strip().rstrip(';'))
    adj = defaultdict(set); idx = {w: i for i, w in enumerate(wit)}; nxt = [len(wit)]
    pos = [0]
    def node():
        if s[pos[0]] == '(':
            me = nxt[0]; nxt[0] += 1; pos[0] += 1
            while True:
                c = node(); adj[me].add(c); adj[c].add(me)
                if s[pos[0]] == ',': pos[0] += 1; continue
                if s[pos[0]] == ')': pos[0] += 1; break
            while pos[0] < len(s) and s[pos[0]] not in ',)': pos[0] += 1  # internal label
            return me
        m = re.match(r"[^,()]+", s[pos[0]:]); name = m.group(0).strip().strip("'"); pos[0] += len(m.group(0))
        if name not in idx:
            me = nxt[0]; nxt[0] += 1; return me  # hypothetical node label
        return idx[name]
    node()
    return adj


def parse_edges(path, wit):
    idx = {w: i for i, w in enumerate(wit)}; names = dict(idx); adj = defaultdict(set)
    def nid(x):
        x = x.strip().strip('"')
        if x not in names: names[x] = len(names)
        return names[x]
    for l in open(path):
        p = l.split()
        if len(p) >= 2:
            a, b = nid(p[0]), nid(p[1]); adj[a].add(b); adj[b].add(a)
    return adj


def clusters(adj, n, ref):
    """For each edge, the witness set on the side not containing ref."""
    out = set()
    for u in list(adj):
        for v in adj[u]:
            seen = {u}; stack = [v]; side = set()
            while stack:
                x = stack.pop()
                if x in seen: continue
                seen.add(x)
                if x < n: side.add(x)
                stack += [y for y in adj[x] if y not in seen]
            if ref not in side and 2 <= len(side) <= n - 2:
                out.add(frozenset(side))
    return out


def score(C, wit, label):
    L = {g: frozenset(wit.index(w) for w in ms) for g, ms in GROUPS.items()}
    rec = sum(s in C for s in L.values())
    prec = sum(any(c == s or c < s for s in L.values()) for c in C)
    jac = {g: round(max((len(c & s) / len(c | s) for c in C), default=0), 3) for g, s in L.items()}
    print(f'{label:12s} families {len(C):3d}  exact {rec}/12  nested-or-equal {prec}/{len(C)}  '
          + ' '.join(f'{g}={"✓" if L[g] in C else jac[g]}' for g in sorted(L, key=lambda g: -len(L[g]))))
    return dict(method=label, families=len(C), recall=rec, precision=prec, jaccard=jac,
                exact={g: L[g] in C for g in L})


if __name__ == '__main__':
    mode, num = sys.argv[1], sys.argv[2]
    js = sys.argv[sys.argv.index('--json') + 1] if '--json' in sys.argv else None
    wit, rows = load_num(num); n = len(wit); ref = wit.index(REF)
    if mode == 'nj': adj, lab = nj_tree(wit, rows), 'NJ'
    elif mode == 'newick': adj, lab = parse_newick(open(sys.argv[3]).read(), wit), os.path.basename(sys.argv[3])
    else: adj, lab = parse_edges(sys.argv[3], wit), os.path.basename(sys.argv[3])
    res = score(clusters(adj, n, ref), wit, lab)
    if js: json.dump(res, open(js, 'w'), indent=1)
