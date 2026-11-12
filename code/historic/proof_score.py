"""Score an extraction against the page-image gold standard (proof/gold_groups.json).
For each page, parsed groups (>=2 sigla, order-insensitive) are matched to printed groups as multisets.
Usage: proof_score.py GOLD.json --units-dir DIR [--mode old|leaf]
  old : units.json 'line' compared with the page's old_lines range (Internet Archive OCR)
  leaf: units.json 'leaf' field equals the page's leaf (re-OCR, one file per leaf)
"""
import json, argparse, os, collections

ap = argparse.ArgumentParser(); ap.add_argument('gold'); ap.add_argument('--units-dir', required=True)
ap.add_argument('--mode', default='old'); a = ap.parse_args()
G = json.load(open(a.gold)); T = collections.Counter()
for ed, pages in G.items():
    if ed.startswith('_'): continue
    f = os.path.join(a.units_dir, ed, 'units.json')
    if not os.path.exists(f): print(ed, 'no units'); continue
    U = json.load(open(f))
    for pg, d in pages.items():
        gold = collections.Counter(''.join(sorted(g)) for g in d['groups'].split())
        if a.mode == 'old':
            lo, hi = d['old_lines']; sel = [u for u in U if lo <= u['line'] < hi]
        else:
            sel = [u for u in U if u.get('leaf') == d['leaf']]
        got = collections.Counter(''.join(sorted(g)) for u in sel for g in u['groups'][1:] if len(g) >= 2)
        right = sum((gold & got).values()); wrong = sum((got - gold).values()); printed = sum(gold.values())
        T['printed'] += printed; T['right'] += right; T['wrong'] += wrong
        print(f'{ed:16s} p.{pg:>4s}  printed {printed:3d}  right {right:3d}  extra {wrong:3d}')
print(f"TOTAL printed {T['printed']}  right {T['right']} ({100*T['right']/T['printed']:.0f}% recall)  extra {T['wrong']}  "
      f"precision {100*T['right']/max(1,T['right']+T['wrong']):.0f}%")
