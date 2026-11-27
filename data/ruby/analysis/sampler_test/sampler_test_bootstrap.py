"""Evidence-side reference for the sampler test: bootstrap over units of WHICH maximal independent set of the register's
families is the maximum-weight one (ties split equally). B = 20,000 resamplings per instance, fixed seeds.

- Priscus: no unit-level file is kept (readers' reports not in repo); the units are rebuilt from the verified tally
  (notes/RUN_PRISCUS_VERIFIED_20261110.md): 45 B+M+P, 8 M+P, 4 E+B, 3 E+P, 3 E+M+P, + 2 places with a dropped family
  (E+B+P, B+P). Each apparatus place supports exactly one rooted family. Sensitivity: + 112 'E alone' places (N=177).
- Lai: unrooted split supports per unit, as code/historic/lai_v2.py (stemma-mwis) computes them. v1 = units.json
  (814 notes; the supports the Ruby register encodes: DE 43) and v2 = units_v2.json (current: DE 41).
- Lawrence: data/lawrence_bhl4759.num + data/lawrence_native_graph_s3.json, unit -> families exactly as
  code/native/bootstrap_histories.py (copied logic; the original is seed 1, B 20000: E 0.7955).
Output: sampler_test_20261123/bootstrap.json"""
import json, os, itertools as it, collections
import numpy as np
from sampler_test_common import INST, OUT, maximal_sets, lab

B = 20000
STEM = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'lawrence'))
res = {}


def boot(name, U, fams, seed, W_full=None):
    """U: list of sets of family names supported by each unit."""
    g = INST[name]; mx = maximal_sets(g); rng = np.random.default_rng(seed); n = len(U)
    X = np.array([[f in u for f in g['fam']] for u in U], dtype=np.int64)
    full = dict(zip(g['fam'], X.sum(0).tolist()))
    M = np.array([[f in s for f in g['fam']] for s in mx], dtype=np.int64)
    P = np.zeros(len(mx)); diffs = []
    for _ in range(B):
        c = np.bincount(rng.integers(0, n, n), minlength=n); w = c @ X; v = M @ w
        top = np.flatnonzero(v == v.max()); P[top] += 1 / len(top)
    P /= B
    vfull = M @ X.sum(0)
    return dict(units=n, B=B, seed=seed, full_weights=full,
                maximal={lab(g, s): dict(size=len(s), weight=int(vfull[k]), bootstrap_share=float(P[k])) for k, s in enumerate(mx)})


# ---------------- Priscus
tally = {'BMP': 45, 'MP': 8, 'EB': 4, 'EP': 3, 'EMP': 3}
U = [{f} for f, k in tally.items() for _ in range(k)] + [set(), set()]
res['priscus'] = boot('priscus', U, None, 1)
res['priscus_sensitivity_N177'] = boot('priscus', U + [set()] * 112, None, 1)
assert res['priscus']['full_weights'] == INST['priscus']['W']

# ---------------- Lai (unrooted split supports, lai_v2.py logic)
LAI = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'lai_ombre')
ALL = frozenset('ABCDEFG')
key = {'CG': 'CG', 'DE': 'DE', 'AB': 'CDEFG', 'DEF': 'DEF', 'EF': 'EF', 'CEG': 'CEG', 'DF': 'DF'}  # register name -> side without A


def unrooted_props(u):
    s = set()
    for gr in u['groups'][1:]:
        if gr and len(gr) >= 2 and frozenset(gr) != ALL:
            g = frozenset(gr)
            if 2 <= len(g) <= 5: s.add(g if 'A' not in g else ALL - g)
    return s


for tag, fn in (('lai_v1_register', 'units.json'), ('lai_v2_current', 'units_v2.json')):
    units = json.load(open(os.path.join(LAI, fn)))
    U = [{r for r, k in key.items() if frozenset(k) in unrooted_props(u)} for u in units]
    res[tag] = boot('lai', U, None, 20261123)
assert res['lai_v1_register']['full_weights'] == INST['lai']['W'], res['lai_v1_register']['full_weights']

# ---------------- Lawrence (bootstrap_histories.py logic)
L = [l for l in open(os.path.join(STEM, 'lawrence_bhl4759.num')).read().split('\n') if l.strip()]
wit = L[0].split(); rows = [l.split() for l in L[1:]]
G = json.load(open(os.path.join(STEM, 'lawrence_native_graph_s3.json')))
famsets = [frozenset(x.split('+')) for x in G['names']]


def unit_fams(r):
    obs = [i for i in range(len(wit)) if r[i] != '?']; out = set()
    for v in set(r[i] for i in obs):
        A = frozenset(wit[i] for i in obs if r[i] == v)
        if 'T' in A: A = frozenset(wit[i] for i in obs) - A
        if A in famsets: out.add(G['names'][famsets.index(A)])
    return out


U = [unit_fams(r) for r in rows]
res['lawrence'] = boot('lawrence', U, None, 1)
assert res['lawrence']['full_weights'] == INST['lawrence']['W'], res['lawrence']['full_weights']
res['lawrence']['reference_original_script'] = 'data/lawrence_native_graph_s3_bootstrap.json: E 0.7955 (python random, seed 1, B 20000)'

json.dump(res, open(os.path.join(OUT, 'bootstrap.json'), 'w'), indent=1)
for k, v in res.items():
    print(k, 'units', v['units'], v['full_weights'])
    for s, r in sorted(v['maximal'].items(), key=lambda t: -t[1]['weight']):
        print(f"   {s:40s} size {r['size']} w {r['weight']:4d}  boot {r['bootstrap_share']:.4f}")
