"""Write the stable families of each classified text as a rooted Newick tree (families must be laminar: nested or disjoint).
Internal nodes are labelled with the family's support (units) and bootstrap stability, e.g. (B,(V,L)'VL w=211 s=1.00').
usage: families_to_newick.py classification_dir out_dir"""
import json, os, sys

def newick(members, fams, label):
    kids = [f for f in fams if f < members and not any(f < g < members for g in fams)]
    covered = set().union(*kids) if kids else set()
    parts = [newick(k, fams, label) for k in sorted(kids, key=lambda f: sorted(f))] + sorted(members - covered)
    return '(' + ','.join(parts) + ')' + label.get(members, '')

cdir, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
for fn in sorted(os.listdir(cdir)):
    if not fn.endswith('.json'): continue
    d = json.load(open(os.path.join(cdir, fn)))
    if 'stable_families' not in d: continue
    W = frozenset(d['analysis_witnesses'])
    fams = [frozenset(k) for k, w, s in d['stable_families']]
    assert all(a <= b or b <= a or not (a & b) for a in fams for b in fams), fn
    label = {frozenset(k): f"'{k} w={w} s={s:.2f}'" for k, w, s in d['stable_families']}
    fams = [f for f in fams if f != W]
    open(os.path.join(out, fn[:-5] + '.nwk'), 'w').write(newick(W, fams, label) + ';\n')
    print(fn[:-5], newick(W, fams, label))
