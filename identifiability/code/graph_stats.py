"""Shape of real stemmatic conflict graphs and the parameters of real traditions.
For each tradition: candidates = non-trivial reading-vs-rest partial splits of the witnesses legible at each unit
(no imputation), weight = number of units; conflict = conflict on the common legible witnesses (overlap conflict).
Reported: n, L, informative units, candidates, conflict edges, density, degrees, share of weight-1 candidates,
components, chordality, the exact MWIS optimum (size, weight share), certificates on the data with the optimum as
reference (share of crossing false/selected pairs with w(false) < w(selected), and whether the local certificate
holds), singular-reading rates per witness (a proxy for pendant rates p_i), the weakest selected family per unit
(a proxy for internal rates), and the long-branch index max_{i,j} p_i p_j against that weakest support.
Data are read from their working locations; nothing is redistributed. Usage: graph_stats.py OUT.json"""
import sys, os, json, itertools, statistics as st
from collections import Counter
import networkx as nx
from ortools.sat.python import cp_model

HOME = os.path.expanduser('~'); P = os.path.join(HOME, 'MC/projects/priscus_excerpta/data')
SRC = {'Heinrichi': os.path.join(HOME, 'stemma_xfer/heinrichi/heinrichi.num'),
       'Donatus': os.path.join(HOME, 'stemma_xfer/donatus_final.num'),
       'Lawrence (BHL 4759)': os.path.join(P, 'lawrence_bhl4759.num'),
       'Jordanes, Getica': os.path.join(P, 'historic/mommsen_jordanes/jordanes.num'),
       'Iwein (Lachmann 1827)': os.path.join(P, 'historic/iwein_1827/iwein.num')}


def load_num(path):
    L = [l.split() for l in open(path).read().split('\n') if l.strip()]
    return L[0], [r for r in L[1:] if len(r) == len(L[0])]


def load_lai():
    wit = list('ABCDEFG'); rows = []
    for u in json.load(open(os.path.join(P, 'historic/lai_ombre/units.json'))):
        r = ['0'] * 7
        for k, g in enumerate(u['groups']):
            for c in g:
                if c in wit: r[wit.index(c)] = str(k)
        rows.append(r)
    return wit, rows


def candidates(wit, rows):
    n = len(wit); w = Counter(); single = Counter(); informative = 0
    for r in rows:
        O = [i for i in range(n) if r[i] != '?']
        props = set()
        for v in set(r[i] for i in O):
            A = frozenset(i for i in O if r[i] == v)
            if len(A) == 1: single[next(iter(A))] += 1
            if 2 <= len(A) <= len(O) - 2:
                B = frozenset(O) - A
                props.add((min(A, B, key=lambda s: min(s)), frozenset(O)))
        informative += bool(props)
        for s in props: w[s] += 1
    return w, single, informative


def conflict(s, t):
    (A, O), (C, Q) = s, t; W = O & Q
    a, b, c, d = A & W, (O - A) & W, C & W, (Q - C) & W
    return bool(a & c and a & d and b & c and b & d)


def mwis(V, w, E):
    m = cp_model.CpModel(); x = {v: m.NewBoolVar('') for v in V}
    for u, v in E: m.Add(x[u] + x[v] <= 1)
    m.Maximize(sum(w[v] * x[v] for v in V))
    s = cp_model.CpSolver(); s.parameters.num_workers = 8; s.parameters.max_time_in_seconds = 600
    st_ = s.Solve(m)
    return [v for v in V if s.Value(x[v])], s.StatusName(st_)


out = {}
data = {k: load_num(v) for k, v in SRC.items() if os.path.exists(v)}
data['Lai de l\'Ombre (Bedier 1913)'] = load_lai()
for name, (wit, rows) in data.items():
    n, L = len(wit), len(rows)
    w, single, informative = candidates(wit, rows)
    V = list(w)
    E = [(u, v) for u, v in itertools.combinations(V, 2) if conflict(u, v)]
    G = nx.Graph(); G.add_nodes_from(V); G.add_edges_from(E)
    deg = [d for _, d in G.degree()]
    S, status = mwis(V, w, E); Sset = set(S)
    cross = [(w[t], w[s]) for t in V if t not in Sset for s in S if G.has_edge(t, s)]
    p = {wit[i]: single[i] / max(1, sum(r[i] != '?' for r in rows)) for i in range(n)}
    top = sorted(p.values(), reverse=True)
    rmin = min(w[s] for s in S) / L if S else None
    out[name] = dict(n=n, L=L, informative_units=informative, candidates=len(V), conflict_edges=len(E),
                     density=round(2 * len(E) / (len(V) * (len(V) - 1)), 4) if len(V) > 1 else 0,
                     degree_median=st.median(deg) if deg else 0, degree_max=max(deg) if deg else 0,
                     weight1_share=round(sum(1 for v in V if w[v] == 1) / len(V), 3),
                     weight_max=max(w.values()), isolated=sum(1 for d in deg if d == 0),
                     components=nx.number_connected_components(G), chordal=nx.is_chordal(G),
                     optimum_size=len(S), optimum_status=status, n_minus_3=n - 3,
                     optimum_weight_share=round(sum(w[s] for s in S) / sum(w.values()), 3),
                     crossing_pairs=len(cross), crossing_pairs_separated=round(sum(a < b for a, b in cross) / max(1, len(cross)), 4),
                     local_certificate=all(a < b for a, b in cross),
                     global_certificate=(min(w[s] for s in S) > max((w[v] for v in V if v not in Sset), default=0)) if S else None,
                     singular_rate_median=round(st.median(p.values()), 4), singular_rate_top2=[round(x, 4) for x in top[:2]],
                     weakest_selected_per_unit=round(rmin, 5) if rmin else None,
                     long_branch_index=round(top[0] * top[1], 5))
    # coincidence factor: for witness pairs that are not a selected family, the rate at which exactly the two
    # share a reading against at least two others (units where both are legible), over the product of their
    # singular-reading rates. Under the irreversible model with a single alternative reading it is about 1.
    sel_pairs = {frozenset(A) for A, O in S if len(A) == 2} | {frozenset(O - A) for A, O in S if len(O - A) == 2}
    kap = []
    for i, j in itertools.combinations(range(n), 2):
        if frozenset({i, j}) in sel_pairs: continue
        both = [r for r in rows if r[i] != '?' and r[j] != '?']
        if not both or p[wit[i]] * p[wit[j]] == 0: continue
        a = sum(1 for r in both if r[i] == r[j] and sum(1 for k in range(n) if r[k] == r[i]) == 2
                and sum(1 for k in range(n) if r[k] not in ('?', r[i])) >= 2) / len(both)
        kap.append((p[wit[i]] * p[wit[j]], a, wit[i], wit[j]))
    kap.sort(reverse=True)
    out[name].update(complete_units=sum(1 for r in rows if '?' not in r),
                     kappa_top_pair=[kap[0][2], kap[0][3]] if kap else None,
                     kappa_top_pair_product=round(kap[0][0], 5) if kap else None,
                     kappa_top_pair_agreement=round(kap[0][1], 5) if kap else None,
                     kappa_top=round(kap[0][1] / kap[0][0], 3) if kap else None,
                     kappa_median=round(st.median(a / pp for pp, a, _, _ in kap), 3) if kap else None,
                     kappa_pairs=len(kap))
    print(name, json.dumps(out[name]))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
