"""Check C21: the metrics of bedier_on_birthdeath.py on hand-made trees, and the claimed orderings in its outputs.
Exits non-zero on any failure. Usage: bedier_birthdeath_check.py RESULT.json [RESULT.json ...]"""
import sys, json, itertools, importlib.util, types
fails = []
def check(ok, msg):
    print(('PASS ' if ok else 'FAIL ') + msg)
    if not ok: fails.append(msg)

# load metrics/prune without running the simulation
src = open(__file__.replace('bedier_birthdeath_check.py', 'bedier_on_birthdeath.py')).read()
mod = types.ModuleType('m'); exec(src.split('results = {}')[0].replace("OUT = sys.argv[1]", "OUT = None").replace("int(sys.argv[2]) if len(sys.argv) > 2 else 300", "0"), mod.__dict__)
F = frozenset
m = mod.metrics(set(), 4)                                    # star on 4 leaves
check(m['i3'] == 0 and not m['bifid_root'] and m['internal'] == 1, 'star: i3 0, root degree 4')
m = mod.metrics({F({0, 1}), F({0, 1, 2})}, 4)                # caterpillar ((0,1),2),3
check(m['i3'] == 1 and m['bifid_root'] and m['d2'] == 3, 'caterpillar: i3 1, all nodes binary')
m = mod.metrics({F({0, 1})}, 4)                              # (0,1),2,3: root degree 3
check(abs(m['i3'] - 0.5) < 1e-12 and m['d3'] == 1 and m['d2'] == 1, '(01)23: i3 1/2, root degree 3')
# prune: root 0 lost with children 1 (alive), 2 (lost, children 3,4 alive), 5 (lost, one child 6 alive)
parent = [-1, 0, 0, 2, 2, 0, 5]; alive = [False, True, False, True, True, False, True]
surv, cl = mod.prune(parent, alive)
check(surv == [1, 3, 4, 6] and cl == {F({1, 2})}, 'prune: lost unary node contracted, lost leaf-free nodes kept as clusters')
# alive internal node: 0 alive with child 1 alive with children 2,3 alive
surv, cl = mod.prune([-1, 0, 1, 1], [True] * 4)
check(cl == {F({1, 2, 3})}, 'prune: surviving exemplar hangs as a leaf at its own node')
check(mod.binom_tail(5, 5) == 1 / 32 and mod.binom_tail(0, 3) == 1, 'binomial tail')

for path in sys.argv[1:]:
    R = json.load(open(path)); tag = path.split('/')[-1]
    for k, r in R.items():
        if not k.startswith('alpha'): continue
        t, raw, riv = r['true'], r['raw'], r['rival']
        check(0.55 <= t['bifid_root'] <= 0.75 and 0.78 <= t['i3'] <= 0.88, f'{tag} {k}: truth in the range of the birth-death null')
        check(raw['bifid_root'] > t['bifid_root'] + 0.1 and raw['i3'] > t['i3'] + 0.05 and raw['p_deg3'] < t['p_deg3'],
              f'{tag} {k}: raw MWIS over-resolves (bifid, i3 up, degree 3 down)')
        rec = r['recovery']; n = t['traditions']
        check(rec['rival']['false_added'] < 0.05 * n, f'{tag} {k}: rival test adds fewer than one false cluster per 20 traditions')
        check(rec['rival']['true_found'] < 0.9 * rec['rival']['true_total'], f'{tag} {k}: rival test loses more than 10% of true clusters')
        check(riv['i3'] < t['i3'] < raw['i3'], f'{tag} {k}: rival under-resolves, raw over-resolves')
        # editors (Camps et al. Table 1, 10% bootstrap CIs): bifid 0.65-0.86, i3 0.89-0.96, P(deg2) 0.76-0.87
        hit = [a for a in ('k>=2', 'k>=3', 'k>=5') if 0.65 <= r[a]['bifid_root'] <= 0.86 and 0.89 <= r[a]['i3'] <= 0.96]
        check(bool(hit), f'{tag} {k}: some support threshold puts bifid and i3 inside the editors\' intervals {hit}')
        check(all(r[a]['p_deg2'] >= 0.85 for a in hit), f'{tag} {k}: the matching arms have P(deg2) >= 0.85, top of the editors\' interval or above {[r[a]["p_deg2"] for a in hit]}')
print('ALL PASS' if not fails else f'{len(fails)} FAILED'); sys.exit(1 if fails else 0)
