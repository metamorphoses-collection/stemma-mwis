"""Compare two independent reader transcriptions (A, B) page by page.
Writes the detailed disagreement list to <out.tsv> (keep it local: it quotes lemmata)
and prints aggregate proof counts only.
Usage: proof_diff.py <dirA> <dirB> <out.tsv>"""
import sys, json, glob, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_matrix import rejoin
from collections import Counter


def key(u):
    lem = re.sub(r'\W+', ' ', u.get('lemma', '').lower()).strip()
    return (str(u.get('chapter')), str(u.get('sentence')), lem)


def sig(u):
    out = []
    for r in u.get('readings', []):
        t = re.sub(r'\W+', ' ', r.get('text', '').lower()).strip()
        out.append((t, tuple(sorted(rejoin(r.get('wits', [])))), r.get('type')))
    return tuple(sorted(out))


A, B, out = sys.argv[1:4]
c = Counter(); rows = []
for pa in sorted(glob.glob(os.path.join(A, 'p*.json'))):
    pb = os.path.join(B, os.path.basename(pa))
    if not os.path.exists(pb):
        c['page_missing_in_B'] += 1; continue
    a, b = json.load(open(pa)), json.load(open(pb))
    ua = {}; ub = {}
    for u in a['units']: ua.setdefault(key(u), []).append(u)
    for u in b['units']: ub.setdefault(key(u), []).append(u)
    c['units_A'] += len(a['units']); c['units_B'] += len(b['units'])
    for k in set(ua) | set(ub):
        if k not in ub:
            c['only_A'] += len(ua[k]); rows += [(a['page'], 'only_A', *k, json.dumps(sig(u), ensure_ascii=False)) for u in ua[k]]
        elif k not in ua:
            c['only_B'] += len(ub[k]); rows += [(a['page'], 'only_B', *k, json.dumps(sig(u), ensure_ascii=False)) for u in ub[k]]
        else:
            sa = sorted(sig(u) for u in ua[k]); sb = sorted(sig(u) for u in ub[k])
            if sa == sb: c['agree'] += len(ua[k])
            else:
                wa = [[w for _, ws, _ in s for w in ws] for s in sa]; wb = [[w for _, ws, _ in s for w in ws] for s in sb]
                kind = 'type_or_text' if wa == wb else 'witnesses'
                c['differ_' + kind] += 1
                rows.append((a['page'], 'differ_' + kind, *k, json.dumps([sa, sb], ensure_ascii=False)))
    ca = sorted((str(x.get('chapter')), str(x.get('sentence')), tuple(sorted(x.get('wits', [])))) for x in a.get('coverage', []))
    cb = sorted((str(x.get('chapter')), str(x.get('sentence')), tuple(sorted(x.get('wits', [])))) for x in b.get('coverage', []))
    c['coverage_agree' if ca == cb else 'coverage_differ'] += 1
    if ca != cb: rows.append((a['page'], 'coverage', '', '', '', json.dumps([ca, cb])))
with open(out, 'w') as f:
    for r in sorted(rows, key=lambda r: (r[0], r[1])):
        f.write('\t'.join(map(str, r)) + '\n')
print(dict(c))
