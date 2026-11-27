"""Scoring for the series test of the editors' top divisions. Implements the rules fixed before any text was classified.

Input: units.json from mgh_rooted_extract.py, and the editor's top division (from data/historic/bedier/tops/*.json).
Steps, as pre-registered:
  1. analysis witnesses = the editor's direct stemma witnesses cited in >= 30 extracted departures (any group size);
  2. non-analysis witnesses are removed from every group; groups equal to all analysis witnesses are dropped;
  3. inclusion: >= 4 analysis witnesses and >= 200 units with a shared departure (>= 2 analysis witnesses);
  4. exact rooted MWIS (min weight 2) and B = 200 bootstrap replicates over units, seed 1 (rooted_mwis.py functions);
     stable family = in the exact optimum with stability >= 0.90;
  5. each editor top branch with >= 2 analysis witnesses: supported (a stable family equals it), contradicted (a stable
     family of weight >= 5 crosses it), else silent; text: supported / contradicted / partly silent;
  6. k_D = maximal stable families + analysis witnesses in no stable family; bootstrap distribution of k_D from each
     replicate's optimum (all its families, not only stable ones).
Usage: bedier_top.py units.json --witnesses ABCD --top 'AB|CD|E' --out result.json [--B 200]
"""
import json, argparse, random, collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rooted_mwis import mwis, weigh

MIN_CITED, MIN_WIT, MIN_UNITS, STABLE, CROSS_W, MINW = 30, 4, 200, 0.90, 5, 2


def crosses(a, b): return bool(a & b) and not (a <= b or b <= a)


def top_blocks(fams, wit):
    fams = [f for f in fams if len(f) >= 2]
    maximal = [f for f in fams if not any(f < g for g in fams)]
    covered = set().union(*maximal) if maximal else set()
    return len(maximal) + len(set(wit) - covered)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('units'); ap.add_argument('--witnesses', required=True)
    ap.add_argument('--top', required=True, help="editor's top branches, e.g. 'HPVLA|SOB|XYZ'")
    ap.add_argument('--out', required=True); ap.add_argument('--B', type=int, default=200)
    ap.add_argument('--screen', action='store_true', help='inclusion counts only; no MWIS (proof pages are scored before classification)')
    a = ap.parse_args()
    U = json.load(open(a.units)); stemma = list(a.witnesses)
    cited = collections.Counter(w for u in U for g in u['groups'][1:] for w in set(g))
    W = [w for w in stemma if cited[w] >= MIN_CITED]; Ws = set(W)
    F = []
    for u in U:
        fs = []
        for g in u['groups'][1:]:
            g = frozenset(x for x in g if x in Ws)
            if len(g) >= 2 and g != Ws: fs.append(g)
        F.append(fs)
    F = [fs for fs in F if fs]
    key = lambda f: ''.join(sorted(f, key=stemma.index))
    res = dict(stemma_witnesses=''.join(stemma), cited={w: cited[w] for w in stemma}, analysis_witnesses=''.join(W),
               units_with_shared_departure=len(F))
    reasons = []
    if len(W) < MIN_WIT: reasons.append('too few witnesses')
    if len(F) < MIN_UNITS: reasons.append('too thin')
    res['included'] = not reasons; res['exclusion'] = reasons
    top = [frozenset(x for x in b if x in Ws) for b in a.top.split('|')]
    top = [b for b in top if b]
    res['editor_top'] = [key(b) for b in top]; res['k_E'] = len(top)
    if reasons or a.screen:
        json.dump(res, open(a.out, 'w'), indent=1); print(json.dumps({k: res[k] for k in ('analysis_witnesses', 'units_with_shared_departure', 'included', 'exclusion')})); return
    Wt = weigh(F, MINW); sel, obj, status = mwis(Wt)
    rng = random.Random(1); stab = collections.Counter(); kd = collections.Counter()
    for _ in range(a.B):
        Fb = [F[rng.randrange(len(F))] for _ in F]
        sb, _, _ = mwis(weigh(Fb, MINW)); stab.update(sb); kd[top_blocks(sb, W)] += 1
    stable = [f for f in sel if stab[f] / a.B >= STABLE]
    branches = []
    for b in top:
        if len(b) < 2: branches.append(dict(branch=key(b), status='single_witness')); continue
        if b in stable: st, why = 'supported', key(b)
        else:
            cr = [f for f in stable if Wt[f] >= CROSS_W and crosses(f, b)]
            st, why = ('contradicted', [key(f) for f in cr]) if cr else ('silent', None)
        branches.append(dict(branch=key(b), status=st, by=why))
    cl = [x['status'] for x in branches if x['status'] != 'single_witness']
    text = 'contradicted' if 'contradicted' in cl else ('supported' if cl and all(s == 'supported' for s in cl) else 'partly_silent')
    res.update(status=status, weight=obj, branches=branches, classification=text,
               k_D=top_blocks(stable, W), k_D_bootstrap={str(k): v / a.B for k, v in sorted(kd.items())},
               stable_families=[[key(f), Wt[f], stab[f] / a.B] for f in sorted(stable, key=lambda f: -Wt[f])],
               mwis=[[key(f), Wt[f], stab[f] / a.B] for f in sorted(sel, key=lambda f: -Wt[f])])
    json.dump(res, open(a.out, 'w'), indent=1)
    print(f"{text}; k_E={len(top)} k_D={res['k_D']}; branches: {branches}")


if __name__ == '__main__':
    main()
