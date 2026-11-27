"""Shared definitions for the check_* scripts. Read-only on the run data.
Graphs are copied from the frozen scorers (score_ruby_mis*.py); register files are the *_mis_prep.json used to score."""
import json, os, itertools as it
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'checks_20261123')
RES = os.path.join(HERE, 'results')

GRAPHS = {
    'priscus': dict(
        fam=['BMP', 'MP', 'EB', 'EP', 'EMP'],
        W={'BMP': 45, 'MP': 8, 'EB': 4, 'EP': 3, 'EMP': 3},
        conflicts=[('BMP', 'EB'), ('BMP', 'EP'), ('BMP', 'EMP'), ('MP', 'EP'), ('EB', 'EP'), ('EB', 'EMP')],
        prep='ruby_mis_prep.json', files=[f'ruby_diag_job_{j}.json' for j in (728, 729, 730, 731)],
        hw_reported=0.62, emu_noisy_reported=0.51),
    'lai': dict(
        fam=['CG', 'DE', 'AB', 'DEF', 'EF', 'CEG', 'DF'],
        W={'CG': 176, 'DE': 43, 'AB': 38, 'DEF': 24, 'EF': 23, 'CEG': 13, 'DF': 13},
        conflicts=[('DE', 'EF'), ('DE', 'CEG'), ('DE', 'DF'), ('DEF', 'CEG'), ('EF', 'CEG'), ('EF', 'DF')],
        prep='lai_mis_prep.json', files=['ruby_diag_job_733.json'],
        hw_reported=0.33, emu_noisy_reported=0.21),
    'lawrence': dict(
        fam=['B+F+Mi+P+R', 'B+R', 'P+R', 'B+F+P+R', 'B+F+R'],
        W={'B+F+Mi+P+R': 11, 'B+R': 4, 'P+R': 4, 'B+F+P+R': 3, 'B+F+R': 3},
        conflicts=[('P+R', 'B+R'), ('P+R', 'B+F+R')],
        prep='lawrence_mis_prep.json', files=['ruby_diag_job_732.json'],
        hw_reported=0.31, emu_noisy_reported=0.14),
}
LAW_E = frozenset({'B+F+Mi+P+R', 'B+R', 'B+F+P+R', 'B+F+R'})
LAW_R = frozenset({'B+F+Mi+P+R', 'P+R', 'B+F+P+R'})


def conflict_fn(g):
    C = {frozenset(e) for e in g['conflicts']}
    return lambda a, b: frozenset((a, b)) in C


def independent_sets(g):
    conf = conflict_fn(g); F = g['fam']
    return [frozenset(s) for k in range(len(F) + 1) for s in it.combinations(F, k)
            if all(not conf(a, b) for a, b in it.combinations(s, 2))]


def maximum_sets(g):
    ind = independent_sets(g); m = max(map(len, ind))
    return [s for s in ind if len(s) == m]


def maximal_sets(g):
    ind = independent_sets(g)
    return [s for s in ind if not any(s < t for t in ind)]


def load_shots(name):
    """Return list of (per-copy outcome tuple, n) for every distinct bitstring, plus file of origin.
    Each bitstring = one shot; the per-copy outcome is the frozenset of families measured in |r>."""
    g = GRAPHS[name]; prep = json.load(open(os.path.join(HERE, g['prep'])))
    order = list(prep['register_um'])
    copies = sorted({q.split('_')[0] for q in order}, key=lambda c: int(c[1:]))
    out = []
    for f in g['files']:
        res = json.load(open(os.path.join(RES, f)))
        assert res['qubit_order'] == order, f'{f}: qubit order differs from the prep register'
        idx = {q: i for i, q in enumerate(order)}
        for bs, n in res['results'].items():
            outs = tuple(frozenset(x for x in g['fam'] if bs[idx[f'{c}_{x}']] == '1') for c in copies)
            out.append((f, outs, int(n)))
    return copies, out


def lab(g, s):
    return '{' + ','.join(x for x in g['fam'] if x in s) + '}'


def dump(obj, name):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name)
    json.dump(obj, open(p, 'w'), indent=1, default=float)
    return p
