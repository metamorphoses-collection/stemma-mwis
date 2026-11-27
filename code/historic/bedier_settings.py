"""Build per-text extraction settings for the Bédier test from the frozen records (tops/*.json), by the rules fixed
before extraction. Mechanical: no choice is made by looking at any result.
Output: data/historic/bedier/settings.json
"""
import json, glob, os, re, random

B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'bedier')
VOLCODE = {'AA 6.1': 61, 'AA 6.2': 62}


def volcode(v):
    if v in VOLCODE: return VOLCODE[v]
    return int(re.search(r'\d+', v).group())


def pages_in(lines, lo, hi):
    """Printed page numbers in OCR lines lo..hi: lone numbers or numbers at either end of a running head, kept only as
    the longest run that increases by 1-4 from one hit to the next (drops years, verse and chapter numbers)."""
    P = []
    for i in range(lo - 1, min(hi, len(lines))):
        t = lines[i].strip()
        m = re.fullmatch(r'(\d{1,3})', t) or re.match(r'^(\d{1,3})\s+[A-Z][A-Z .,]{5,}$', t) or re.match(r'^[A-Z][A-Z .,]{5,}\s+(\d{1,3})$', t)
        if m: P.append(int(m.group(1)))
    best = []                                   # longest chain p1 < p2 < ... with gaps 1..4, in OCR order
    chains = []
    for p in P:
        ext = [c + [p] for c in chains if 1 <= p - c[-1] <= 4]
        chains = [c for c in chains if not (1 <= p - c[-1] <= 4)] + (max(ext, key=len) and [max(ext, key=len)] if ext else [[p]])
        chains = sorted(chains, key=len, reverse=True)[:50]
    if chains: best = max(chains, key=len)
    if not best: return []
    return list(range(best[0], best[-1] + 1))   # every page between first and last found


out = []
for f in sorted(glob.glob(os.path.join(B, 'tops', '*.json'))):
    d = json.load(open(f)); lines = open(os.path.join(B, 'ocr', d['ia_id'] + '.txt'), encoding='utf-8', errors='replace').readlines()
    for idx, t in enumerate(d['texts']):
        if t['classification_A3'] not in ('explicit', 'implied', 'undecided'): continue
        direct = [s['siglum'] for s in t['sigla'] if s.get('role') == 'direct']
        caps = [s for s in direct if re.fullmatch(r'[A-Z]', s)]
        caps = list(dict.fromkeys(caps))
        digits = 'drop' if any(re.fullmatch(r'[A-Z]\d', s) for s in direct) else 'hand'
        fixes = ','.join(x for x, c in (('0:O', 'O'), ('8:S', 'S')) if c in caps)
        top = ['' .join(s for s in b if s in caps) for b in t.get('top_branches', [])]
        top = [b for b in top if b]
        lo, hi = t['text_lines']
        pg = pages_in(lines, lo, hi)
        rng = random.Random(volcode(d['volume']) * 1000 + idx)
        proof = rng.choice(pg) if pg else None
        out.append(dict(volume=d['volume'], ia_id=d['ia_id'], text_id=t['text_id'], index=idx, kind=t['classification_A3'],
                        pilot='pilot_text' in t.get('flags', []), lines=[lo, hi], sigla=''.join(caps),
                        other_direct=[s for s in direct if s not in caps], digits=digits, ocr=fixes,
                        top='|'.join(top), k_E_recorded=len(t.get('top_branches', [])), n_pages=len(pg), proof_page=proof))
json.dump(out, open(os.path.join(B, 'settings.json'), 'w'), ensure_ascii=False, indent=1)
for s in out:
    print(f"{s['text_id'][:30]:30s} {s['kind'][:4]} {'P' if s['pilot'] else ' '} sig={s['sigla']:14s} top={s['top']:28s} dig={s['digits']} pages={s['n_pages']:3d} proof={s['proof_page']}")
