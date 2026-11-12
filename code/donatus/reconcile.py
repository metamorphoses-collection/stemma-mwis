"""Build the proofed transcription: reader A as base (fuller hand / '(om. X)' encoding) plus the
fixes decided on the page images after the A/B diff. Also reports matrix-level A/B agreement.
Fixes live in a LOCAL json (keyed by page, chapter, sentence and lemma prefix); they only mark
witnesses as unknown ('?'). No lemma text is kept in this repo.
Usage: reconcile.py <dirA> <dirB> <outdir>"""
import sys, os, json, glob, re, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_matrix import build, load, WIT

# (page, chapter, sentence, lemma prefix, witnesses whose state is unknown, reason)
FIXES_FILE = os.environ.get('PROOF_FIXES', 'proof_fixes.json')  # local: keyed on lemmata (copyright)
FIXES = [(f['page'], f['chapter'], f['sentence'], f['lemma_prefix'], f['unknown'], f['reason']) for f in json.load(open(FIXES_FILE))]

A, B, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
applied = 0
for f in sorted(glob.glob(os.path.join(A, 'p*.json'))):
    d = json.load(open(f))
    for (pg, ch, se, pre, wits, why) in FIXES:
        if d['page'] != pg: continue
        for u in d['units']:
            if str(u.get('chapter')) == ch and (se is None or str(u.get('sentence')) == se) \
                    and u.get('lemma', '').lower().startswith(pre):
                u['readings'].append({'text': '[proof: unknown]', 'wits': wits, 'type': 'other'})
                applied += 1
                break
    json.dump(d, open(os.path.join(out, os.path.basename(f)), 'w'), ensure_ascii=False)
print(f'fixes applied {applied}/{len(FIXES)}')


# matrix-level agreement A vs B on units both readers have (same chapter, sentence, lemma, order)
def keyed(d):
    rows, man, *_ = build(load(d))
    seen = {}; out = {}
    for r, m in zip(rows, man):
        k = (m[1], m[2]); seen[k] = seen.get(k, 0) + 1
        out[(k, seen[k])] = r
    return out
ka, kb = keyed(A), keyed(B)
common = set(ka) & set(kb)
cells = agree = 0
for k in common:
    for x, y in zip(ka[k], kb[k]):
        cells += 1; agree += (x == y) or ('?' in (x, y) and False)
print(f'rows A {len(ka)}  B {len(kb)}  aligned {len(common)}  cell agreement {agree}/{cells} = {agree / cells:.4f}')
