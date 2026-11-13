"""Score extraction quality on the pre-drawn proof page of each included text (PREREG rule 5: precision >= 85%,
recall >= 60%). Gold = groups transcribed from the dMGH page image as printed; the counting rule is applied here, in
code, identically for every text:
  - a sigla string is split into sigla; hands: X or X1 = first hand (kept as X), X2+ = corrector (dropped), except for
    texts with digits='drop' (A4), where letter+digit sigla are other manuscripts and are dropped;
  - only sigla of the text's analysis form (settings 'sigla': single capitals) are kept;
  - a reading's group counts if it keeps >= 2 sigla (lemma sigla are the editor's text: not counted, as in the pilot).
Parsed groups = groups of >= 2 in the units whose OCR line falls on that page (between the page-number line of the page
and of the next page in the text's range). Matching is order-insensitive, as multisets.
Usage: bedier_proof.py GOLD.json [GOLD2.json ...] --leaves leaves.json --out proof_scores.json
"""
import json, re, os, argparse, collections

B = os.path.expanduser('~/MC/projects/priscus_excerpta/data/historic/bedier')


def split_sigla(s):
    """Tokens of a transcribed sigla string, whatever the reader's spacing: runs like 'VF' or 'MbTaCFP' are split into
    sigla, each a capital with an optional lower-case or numeric superscript ('Mb', 'P2'); other words ('pr.', 'm.') drop."""
    s = s.replace('¹', '1').replace('²', '2').replace('³', '3')
    toks = []
    for w in s.split():
        if re.fullmatch(r'(?:[A-Z][a-z]?\d?)+', w) and not re.fullmatch(r'[A-Z][a-z]{2,}\.?', w):
            toks += re.findall(r'[A-Z][a-z]?\d?', w)
        else:
            toks.append(w)
    return toks


def gold_groups(entry, sig, digits, lower=''):
    G = []
    for r in entry.get('readings', []):
        keep = []
        for tok in split_sigla(r.get('sigla', '')):
            ml = re.fullmatch(r'([A-Z])([a-z])(\d?)', tok)              # A10: same lower-case rule as the parser
            if ml and ml.group(2) in lower: tok = ml.group(1) + ml.group(3)
            m = re.fullmatch(r'([A-Z])(\d*)', tok)
            if not m: continue
            letter, d = m.groups()
            if d and digits == 'drop': continue
            if d and d != '1': continue
            if letter in sig: keep.append(letter)
        keep = sorted(set(keep))
        if len(keep) >= 2: G.append(''.join(keep))
    return G


def _pg(t):
    """Page number on a page-head line: a short line whose first or last token is a 1-3 digit number, alone or beside a
    running head (any word of >= 4 letters, however garbled by OCR)."""
    t = t.strip()
    if len(t) > 70: return None
    m = re.fullmatch(r'(\d{1,3})', t)
    if m: return int(m.group(1))
    w = t.split()
    if len(w) >= 2 and re.search(r'[A-Za-z]{4,}', t):
        for tok in (w[0], w[-1]):
            m = re.fullmatch(r'[:.,]?(\d{1,3})[.,]?', tok)
            if m: return int(m.group(1))
    return None


def page_lines(lines, lo, hi, p):
    """OCR window of printed page p: from its page-head line to the next page-head line of p+1..p+3; if none is found,
    the text's median page length (from consecutive page heads found in the text's range)."""
    heads = [(i + 1, n) for i in range(lo - 1, min(hi, len(lines))) for n in [_pg(lines[i])] if n is not None]
    gaps = sorted(b[0] - a[0] for a, b in zip(heads, heads[1:]) if 1 <= b[1] - a[1] <= 2 and b[0] > a[0])
    med = gaps[len(gaps) // 2] // max(1, 1) if gaps else 100
    a = next((i for i, n in heads if n == p), None)
    if a is None: return None
    b = next((i for i, n in heads if i > a and p < n <= p + 3), None)
    return a, (b if b is not None and b - a <= 2 * med else a + med)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('gold', nargs='+'); ap.add_argument('--out', required=True); ap.add_argument('--windows', help='A5: content-aligned page windows (bedier_align.py)')
    ap.add_argument('--leaves', help='A6/A7: proof leaves (leaves.json); units read from --runs, matched by leaf')
    ap.add_argument('--runs', default=os.path.join(B, 'runs'))
    a = ap.parse_args()
    st = {s['text_id']: s for s in json.load(open(os.path.join(B, 'settings.json')))}
    res = {}
    WIN = json.load(open(a.windows)) if a.windows else None
    for g in a.gold:
        for pg in json.load(open(g))['pages']:
            t = pg['text_id']; s = st[t]
            gold = collections.Counter(x for e in pg['entries'] for x in gold_groups(e, set(s['sigla']), s['digits'], s.get('lower', '')))
            if a.leaves:
                LV = json.load(open(a.leaves)); pl = LV[t]; rf = os.path.join(a.runs, t, 'units.json')
                if not os.path.exists(rf): res[t] = dict(page=pg['page'], error='not extracted (not located or volume not OCRed)'); continue
                rl = json.load(open(os.path.join(a.runs, t, 'leaves.json')))
                if rl['vol'] != pl['vol'] or not (rl['leaves'][0] <= pl['leaf'] <= rl['leaves'][1]):
                    res[t] = dict(page=pg['page'], error=f"proof leaf {pl['vol']} {pl['leaf']} outside the text's leaves {rl}"); continue
                U = json.load(open(rf))
                got = collections.Counter(''.join(sorted(gr)) for u in U if u['leaf'] == pl['leaf'] for gr in u['groups'][1:] if len(gr) >= 2)
                right = sum((gold & got).values()); extra = sum((got - gold).values()); printed = sum(gold.values())
                rec = right / printed if printed else None; prec = right / (right + extra) if right + extra else None
                if printed == 0: ok = extra == 0
                else: ok = prec is not None and rec >= 0.60 and prec >= 0.85
                res[t] = dict(page=pg['page'], leaf=pl['leaf'], printed=printed, right=right, extra=extra, recall=rec, precision=prec, passes=ok,
                              missed=sorted((gold - got).elements()), extras=sorted((got - gold).elements()))
                print(f"{t[:32]:32s} p.{pg['page']:>4}  printed {printed:3d} right {right:3d} extra {extra:3d}  "
                      f"rec {rec if rec is None else round(rec, 2)}  prec {prec if prec is None else round(prec, 2)}  {'PASS' if ok else 'FAIL'}")
                continue
            lines = open(os.path.join(B, 'ocr', s['ia_id'] + '.txt'), encoding='utf-8', errors='replace').readlines()
            if WIN is not None:
                w = WIN.get(t, {}); rng = (w['start'], w['end']) if w.get('start') and w.get('end') and w['end'] > w['start'] else None
            else:
                rng = page_lines(lines, s['lines'][0], s['lines'][1], pg['page'])
            U = json.load(open(os.path.join(B, 'runs', t, 'units.json')))
            if rng is None:
                res[t] = dict(page=pg['page'], error='page not found in OCR range'); continue
            lo, hi = rng
            got = collections.Counter(''.join(sorted(gr)) for u in U if lo <= u.get('entry_line', u['line']) < hi for gr in u['groups'][1:] if len(gr) >= 2)
            right = sum((gold & got).values()); extra = sum((got - gold).values()); printed = sum(gold.values())
            rec = right / printed if printed else None; prec = right / (right + extra) if right + extra else None
            if printed == 0: ok = extra == 0        # fixed before any gold was seen: an empty page passes only if nothing is parsed
            else: ok = prec is not None and rec >= 0.60 and prec >= 0.85
            res[t] = dict(page=pg['page'], ocr_lines=[lo, hi], printed=printed, right=right, extra=extra,
                          recall=rec, precision=prec, passes=ok, missed=sorted((gold - got).elements()), extras=sorted((got - gold).elements()))
            print(f"{t[:32]:32s} p.{pg['page']:>4}  printed {printed:3d} right {right:3d} extra {extra:3d}  "
                  f"rec {rec if rec is None else round(rec, 2)}  prec {prec if prec is None else round(prec, 2)}  {'PASS' if ok else 'FAIL'}")
    json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
