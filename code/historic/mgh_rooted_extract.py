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


def make(sigla, ocr, digits='hand', lower=''):
    S = set(sigla); LOW = set(lower)
    fix = dict(p.split(':') for p in ocr.split(',') if p) if ocr else {}

    def sig_of(tok):
        """Return (sigla, comma) for a sigla token. Hands: X or X1 = first hand (kept), X2.. = corrector (dropped).
        OCR repairs (v2): hand marks rendered as ! ? ^ ' are read as 1 / 2 / 1 / 1; 'l' read as I where I is a
        siglum, else as the hand mark 1; 'Q' read as G where G is a siglum and Q is not; a doubled letter (TT) is one."""
        t = tok.strip(',.;:|')
        if not t or len(t) > 16: return None, False
        if LOW and len(t) <= 4 and re.fullmatch(r'(?:[a-z]\d?)+', t) and set(re.findall('[a-z]', t)) <= LOW and len(set(t)) == len(re.findall('[a-z]', t)):
            return '', tok.rstrip().endswith(',')                       # v4: standalone lower-case siglum (non-analysis)
        u = ''.join(fix.get(c, c) for c in t)
        u = u.replace('!', '1').replace('¹', '1').replace('²', '2').replace('?', '2').replace('^', '1').replace("'", '1').replace('’', '1')
        if 'I' in S: u = re.sub(r'(?<=[A-Z])l', 'I', u)
        else: u = re.sub(r'(?<=[A-Z])l', '1', u)
        if 'I' not in S: u = re.sub(r'(?<=[A-Z])I(?![A-Z]?[a-z])', '1', u)   # v3: hand mark 1 read as capital I (Vision)
        if 'G' in S and 'Q' not in S: u = u.replace('Q', 'G')
        u = re.sub(r'([A-Z])\1', r'\1', u)
        # v4 (A10): a lower-case letter after a capital is a separate non-analysis siglum if the text's record lists it
        # (dropped, the capital kept), else a superscript (that siglum dropped, as in the scorer)
        u = re.sub(r'([A-Z]\d?)([a-z])', lambda m: m.group(1) if m.group(2) in LOW else '#', u).replace('#', '')
        if not u: return '', tok.rstrip().endswith(',')
        if not re.fullmatch(r'(?:[A-Z]\d?)+', u): return None, False
        # v4 (A10): capitals outside the text's sigla inside a run of >= 2 letters are other witnesses: dropped
        if len(re.findall(r'[A-Z]', u)) >= 2 and any(c in S for c in re.findall(r'[A-Z]', u)):
            u = re.sub(r'[A-Z]\d?', lambda m: m.group(0) if m.group(0)[0] in S else '', u)
        if digits == 'drop' and re.search(r'\d', u): return '', tok.rstrip().endswith(',')   # A4: digit = other manuscript, not a hand
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
    if re.search(r'(^|\s)[a-z]\)\s', t) and n >= 1: return True                   # v4: footnotes keyed a) b) ...
    w = t.split()                                                              # v4: interlinear apparatus (no lemma, no line no.)
    pairs = sum(1 for a_, b_ in zip(w, w[1:]) if re.fullmatch(r"[a-z][\w.'-]*", a_) and b_[-1:] in ',;:' and sig_of(b_)[0])
    if pairs >= 2 or (re.match(r'^om\.\s', t) and n >= 1): return True
    if len(w) == 2 and re.fullmatch(r"[a-z][\w'-]*", w[0]) and sig_of(w[1])[0]: return True
    return in_block and n >= 1


def entries(text, sig_of):
    text = re.sub(r'(\w)-\s+(@@L\d+)\s+(\w\S*)', r'\1\3 \2', text)       # hyphen across a line marker: move the whole token
    text = re.sub(r'(\w)-\s+(\w)', r'\1\2', text)
    text = re.sub(r'\([^()]*\)', ' ', text.replace('||', ' '))
    out = []; cur = None; buf = []

    cur_line = [None]

    def flush():
        nonlocal cur
        if cur and cur['groups']: out.append(cur)
        cur = None

    toks = text.split()
    merged = []
    for t in toks:                                  # "M TCFP" -> "MTCFP" when both are sigla and letters stay distinct
        if t.startswith('@@L'): merged.append(t); continue
        j = next((i for i in range(len(merged) - 1, -1, -1) if not merged[i].startswith('@@L')), None)
        if j is not None and not merged[j].endswith(',') and sig_of(merged[j])[0] and sig_of(t)[0]:
            a_, b_ = sig_of(merged[j])[0], sig_of(t)[0]
            if not set(a_) & set(b_) and len(a_) + len(b_) <= 26:
                merged[j] = merged[j] + t; continue
        merged.append(t)
    toks = merged
    real = [i for i, t in enumerate(toks) if not t.startswith('@@L')]          # neighbours ignore line markers
    nxt = {a: b for a, b in zip(real, real[1:])}; prv = {b: a for a, b in zip(real, real[1:])}
    caps = lambda w: bool(re.fullmatch(r'[A-Z][A-Z0-9€]{2,}[.,;:]?', w)) and len(re.findall('[A-Z]', w)) >= 3 and not sig_of(w)[0]
    for k, tok in enumerate(toks):
        if tok.startswith('@@L'):
            cur_line[0] = int(tok[3:]); continue
        # title guard: a sigla-like token beside a capitalised word (inscriptions: "ITEM AD BVNDEM") is text
        if sig_of(tok)[0] and ((k in nxt and caps(toks[nxt[k]])) or (k in prv and caps(toks[prv[k]]))):
            buf.append(tok); continue
        if re.fullmatch(r'\d{1,3}[.,]?', tok) or re.fullmatch(r"[a-z]'?\)", tok):
            flush(); buf = []; continue
        if tok.endswith(']'):
            flush(); cur = dict(lemma=' '.join(buf + [tok[:-1]]), lemma_sig=None, groups=[], eline=cur_line[0]); buf = []; continue
        s, comma = sig_of(tok)
        if s == '':
            buf = []
            if not comma: flush()
            continue
        if s:
            if cur is None: cur = dict(lemma=None, lemma_sig=None, groups=[], eline=cur_line[0])
            if not buf and cur['lemma'] is not None and not cur['groups'] and cur['lemma_sig'] is None:
                cur['lemma_sig'] = s
            elif tok.rstrip().endswith(':') and not cur['groups'] and cur['lemma_sig'] is None:
                cur['lemma'] = ' '.join(buf); cur['lemma_sig'] = s      # v4: "lemma SIGLA: reading SIGLA"
                buf = []; continue
            elif buf:
                used = set(cur['lemma_sig'] or '') | {w for _, g in cur['groups'] for w in g}
                if used & set(s):                                       # v4: a repeated siglum starts a new entry
                    flush(); cur = dict(lemma=None, lemma_sig=None, groups=[], eline=cur_line[0])
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
    ap.add_argument('--lower', default='', help='A10: lower-case letters that are separate sigla in this edition (from the blind record)')
    ap.add_argument('--digits', default='hand', choices=['hand', 'drop'], help="A4: 'drop' where letter+digit sigla are distinct manuscripts")
    ap.add_argument('--strict-blocks', action='store_true', help='v1 behaviour: any non-apparatus line closes a block')
    ap.add_argument('--skip-from', default='', help='regex: start of a region to skip (e.g. a second text on the page)')
    ap.add_argument('--skip-to', default='', help='regex: end of the skipped region (e.g. the next running head)')
    a = ap.parse_args()
    sig_of = make(a.sigla, a.ocr, a.digits, a.lower)
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
                U.append(dict(line=start, entry_line=e.get('eline') or start, leaf=cur_leaf[0], lemma=e['lemma'], lemma_group=0,
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
            block.append(f'@@L{i + 1} ' + t)
        elif t:
            nxt = [L[j].strip() for j in range(i + 1, min(i + 6, hi)) if L[j].strip()][:2]
            if start is not None and any(is_app(x, sig_of, True) for x in nxt) and not a.strict_blocks:
                block.append(f'@@L{i + 1} ' + t)         # one garbled line inside an apparatus block
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
