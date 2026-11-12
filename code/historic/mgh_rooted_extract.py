"""Generic extraction of shared departures from a negative apparatus (MGH style), for the rooted analysis.

Works without page headers: only groups of witnesses that share a reading departing from the edited text are kept
(Lachmann's rule), so absent witnesses never need to be known.  Entry forms handled:
    lemma] reading SIGLA, reading SIGLA        -> every cited group after ']' is a departure,
    lemma] SIGLA, reading SIGLA                   except a first group with no reading (the lemma's own sigla)
    reading SIGLA                                 -> departure
    word om. SIGLA                                -> departure (omission)
An apparatus line is a line containing ']' or starting with a line number and citing at least two sigla tokens.
Parenthesised sub-variants are dropped; tokens with digits (correctors, hands) are ignored; a siglum cited twice
in one entry drops the entry.  Output: units.json in the format read by rooted_mwis.py.

Usage: mgh_rooted_extract.py <djvu.txt> --sigla ABCD... --out DIR [--from N] [--to N] [--name NAME]
       [--ocr '0:O,8:S']   (character fixes applied inside sigla tokens only)
"""
import re, json, argparse, os
from collections import Counter


def make(sigla, ocr):
    S = set(sigla)
    fix = dict(p.split(':') for p in ocr.split(',') if p) if ocr else {}

    def sig_of(tok):
        """Return (sigla, comma) for a sigla token. Hands: X or X1 = first hand (kept), X2.. = corrector (dropped).
        OCR repairs (v2): hand marks rendered as ! ? ^ ' are read as 1 / 2 / 1 / 1; 'l' read as I where I is a
        siglum, else as the hand mark 1; 'Q' read as G where G is a siglum and Q is not; a doubled letter (TT) is one."""
        t = tok.strip(',.;:|')
        if not t or len(t) > 16: return None, False
        u = ''.join(fix.get(c, c) for c in t)
        u = u.replace('!', '1').replace('¹', '1').replace('²', '2').replace('?', '2').replace('^', '1').replace("'", '1').replace('’', '1')
        if 'I' in S: u = re.sub(r'(?<=[A-Z])l', 'I', u)
        else: u = re.sub(r'(?<=[A-Z])l', '1', u)
        if 'G' in S and 'Q' not in S: u = u.replace('Q', 'G')
        u = re.sub(r'([A-Z])\1', r'\1', u)
        if not re.fullmatch(r'(?:[A-Z]\d?)+', u): return None, False
        out = ''.join(m.group(1) for m in re.finditer(r'([A-Z])(\d?)', u) if m.group(2) in ('', '1'))
        letters = re.findall(r'[A-Z]', u)
        if not set(letters) <= S or len(set(letters)) != len(letters): return None, False
        if not out: return '', tok.rstrip().endswith(',')          # only correctors: a sigla token with no group
        return out, tok.rstrip().endswith(',')
    return sig_of


ROMAN = re.compile(r'^(?=[MDCLXVI]{2,}$)M{0,3}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$')


def sig_count(line, sig_of):
    return sum(1 for w in line.split() if sig_of(w)[0] is not None and not ROMAN.match(w.strip(',.;:|')))


def is_app(line, sig_of, in_block=False):
    t = line.strip()
    if not t: return False
    n = sig_count(t, sig_of)
    if n >= 1 and (']' in t or '||' in t or ' om.' in t): return True
    if re.match(r'^\d{1,2}\s', t) and n >= 2: return True
    return in_block and n >= 1


def entries(text, sig_of):
    text = re.sub(r'(\w)-\s+(\w)', r'\1\2', text)
    text = re.sub(r'\([^()]*\)', ' ', text.replace('||', ' '))
    out = []; cur = None; buf = []

    def flush():
        nonlocal cur
        if cur and cur['groups']: out.append(cur)
        cur = None

    toks = text.split()
    merged = []
    for t in toks:                                  # "M TCFP" -> "MTCFP" when both are sigla and letters stay distinct
        if merged and not merged[-1].endswith(',') and sig_of(merged[-1])[0] and sig_of(t)[0]:
            a_, b_ = sig_of(merged[-1])[0], sig_of(t)[0]
            if not set(a_) & set(b_) and len(a_) + len(b_) <= 26:
                merged[-1] = merged[-1] + t; continue
        merged.append(t)
    toks = merged
    caps = lambda w: bool(re.fullmatch(r'[A-Z][A-Z0-9€]{2,}[.,;:]?', w)) and len(re.findall('[A-Z]', w)) >= 3 and not sig_of(w)[0]
    for k, tok in enumerate(toks):
        # title guard: a sigla-like token beside a capitalised word (inscriptions: "ITEM AD BVNDEM") is text
        if sig_of(tok)[0] and ((k + 1 < len(toks) and caps(toks[k + 1])) or (k and caps(toks[k - 1]))):
            buf.append(tok); continue
        if re.fullmatch(r'\d{1,3}[.,]?', tok):
            flush(); buf = []; continue
        if tok.endswith(']'):
            flush(); cur = dict(lemma=' '.join(buf + [tok[:-1]]), lemma_sig=None, groups=[]); buf = []; continue
        s, comma = sig_of(tok)
        if s == '':
            buf = []
            if not comma: flush()
            continue
        if s:
            if cur is None: cur = dict(lemma=None, lemma_sig=None, groups=[])
            if not buf and cur['lemma'] is not None and not cur['groups'] and cur['lemma_sig'] is None:
                cur['lemma_sig'] = s
            elif buf:
                cur['groups'].append((' '.join(buf), s))
            buf = []
            if not comma: flush()
            continue
        buf.append(tok)
    flush()
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('txt'); ap.add_argument('--sigla', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--from', dest='lo', type=int, default=1)
    ap.add_argument('--to', dest='hi', type=int, default=0); ap.add_argument('--ocr', default='')
    ap.add_argument('--name', default='edition')
    ap.add_argument('--strict-blocks', action='store_true', help='v1 behaviour: any non-apparatus line closes a block')
    ap.add_argument('--skip-from', default='', help='regex: start of a region to skip (e.g. a second text on the page)')
    ap.add_argument('--skip-to', default='', help='regex: end of the skipped region (e.g. the next running head)')
    a = ap.parse_args()
    sig_of = make(a.sigla, a.ocr)
    L = open(a.txt, encoding='utf-8', errors='replace').readlines()
    hi = a.hi or len(L)
    U = []; block = []; start = None; cur_leaf = [None]

    def close():
        nonlocal block, start
        if block:
            for e in entries(' '.join(block), sig_of):
                if (e['lemma'] or '').startswith('subscr') or any(r.startswith('subscr') for r, _ in e['groups']):
                    continue                                    # stichometric subscriptions, not readings
                gs = [g for _, g in e['groups']]
                allc = [w for g in gs for w in g] + list(e['lemma_sig'] or '')
                if len(allc) != len(set(allc)): continue
                U.append(dict(line=start, leaf=cur_leaf[0], lemma=e['lemma'], lemma_group=0,
                              groups=[e['lemma_sig'] or ''] + gs, readings=[r for r, _ in e['groups']]))
        block = []; start = None

    skip = False
    for i in range(a.lo - 1, hi):
        t = L[i].strip()
        mleaf = re.match(r'^@@LEAF (\d+)', t)
        if mleaf:
            close(); cur_leaf[0] = int(mleaf.group(1)); skip = False; continue
        if a.skip_from and re.search(a.skip_from, t): close(); skip = True; continue
        if skip:
            if a.skip_to and re.search(a.skip_to, t): skip = False
            continue
        if is_app(t, sig_of, in_block=start is not None):
            if start is None: start = i + 1
            block.append(t)
        elif t:
            nxt = [L[j].strip() for j in range(i + 1, min(i + 6, hi)) if L[j].strip()][:2]
            if start is not None and any(is_app(x, sig_of, True) for x in nxt) and not a.strict_blocks:
                block.append(t)                          # one garbled line inside an apparatus block
            else:
                close()
    close()
    os.makedirs(a.out, exist_ok=True)
    json.dump(U, open(os.path.join(a.out, 'units.json'), 'w'), ensure_ascii=False, indent=0)
    fam = Counter(g for u in U for g in u['groups'][1:] if len(g) >= 2)
    print(f'{a.name}: {len(U)} entries, {sum(fam.values())} shared departures in {len(fam)} groups')
    print('top:', fam.most_common(20))


if __name__ == '__main__':
    main()
