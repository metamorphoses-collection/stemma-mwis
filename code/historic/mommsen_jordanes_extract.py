"""Deterministic extraction of Mommsen's Getica apparatus (MGH Auct. ant. 5,1, 1882) into a partial collation matrix.

Mommsen's apparatus is a negative apparatus in the form `lemma] reading SIGLA, reading SIGLA` or
`reading SIGLA`, page by page, under a header listing the witnesses present on that page,
e.g. `H P V L A O B X (Y Z).`  Parenthesised witnesses (Y, Z) and those with inc./fin. notes are
reported only selectively: they are treated as positive-only (state known only where cited).

Each entry becomes one variation unit: each cited group of sigla is one state; witnesses present
on the page (non-parenthesised) and not cited get the lemma state, unless the lemma carries its own
explicit sigla, in which case uncited witnesses are unknown ('?').  No LLM; every unit is traceable
to the OCR line it came from.

Usage: mommsen_jordanes_extract.py <djvu.txt> --out DIR [--from LINE] [--to LINE]
"""
import re, json, argparse, os
from collections import Counter

SIG = 'HPVLAOBXYZS'
ORDER = list('HPVLASOBXYZ')
FULL = set('HPVLAOBX')
HEADER = re.compile(r'^\s*H\s?P\s?V\s?L')
SIGTOK = re.compile(r'^[HPVLAOBXYZS0]+$')


def fix_header(s):
    return s.replace('{', '(').replace('}', ')').replace('[', '(').replace(']', ')').replace('7', 'Y').replace('8', 'S').replace('0', 'O')


def parse_header(line):
    """Return (core, partial) sigla sets, or None if the line is an apparatus line starting with sigla."""
    s = fix_header(line.strip())
    m = re.match(r'^([HPVLAOBXYZS\s]+(?:\([^)]*\)[HPVLAOBXYZS\s]*)*)', s)
    if not m: return None
    head = m.group(1)
    rest = s[m.end():]
    if rest.startswith(','): return None          # e.g. "HPVLOBXY, bitiniae A": apparatus, not header
    core, partial = set(), set()
    depth = 0; note = ''
    for ch in head:
        if ch == '(': depth += 1; note = ''
        elif ch == ')': depth -= 1
        elif ch in SIG:
            (partial if depth else core).add(ch)
    # witnesses with "(inc ...)" / "(fin ...)" notes right after them are partial
    # a witness followed by a range or inc./fin. note is partial; "X (Y Z)" does not make X partial
    for m2 in re.finditer(r'([HPVLAOBXYZS])\s*\((?![HPVLAOBXYZS][\s)])', s):
        core.discard(m2.group(1)); partial.add(m2.group(1))
    return core, partial


def is_stop(line):
    t = line.strip()
    w = re.match(r'^([A-Z]{5,})\b', t)
    running_head = bool(w and set(w.group(1)) - set(SIG))      # GETICA, IORDANIS…; not a sigla run like VABXY
    return bool(re.match(r'^\d+\)\s', t) or re.match(r'^\d+\s+[A-Z]{4,}', t) or running_head
                or re.match(r'^[IVXL]+\s+\d+\s', t))


def blocks(lines, lo, hi):
    """Yield (header_line_no, core, partial, [(line_no, text)]) apparatus blocks."""
    i = lo
    while i < hi:
        if HEADER.match(lines[i]):
            h = parse_header(lines[i])
            if h:
                body = []; j = i + 1
                while j < hi and not is_stop(lines[j]) and not (HEADER.match(lines[j]) and parse_header(lines[j])):
                    if lines[j].strip(): body.append((j + 1, lines[j].rstrip('\n')))
                    j += 1
                yield i + 1, h[0], h[1], body
                i = j; continue
        i += 1


def sig_of(tok):
    t = tok.strip(',.;:')
    t = re.sub(r'[\^]', '', t)
    if t and SIGTOK.match(t) and (t.replace('0', '') != '' or t == '0'):
        t = t.replace('0', 'O')
        if len(set(t)) == len(t): return t, tok.rstrip().endswith(',')
    return None, False


def entries(text):
    text = re.sub(r'(\w)-\s+(\w)', r'\1\2', text)            # rejoin hyphenated line breaks
    text = re.sub(r'\([^()]*\)', ' ', text)                 # drop parenthesised sub-variants
    toks = text.split()
    out = []; cur = None; buf = []

    def flush():
        nonlocal cur
        if cur and cur['groups']: out.append(cur)
        cur = None

    for tok in toks:
        if re.fullmatch(r'\d{1,2}', tok):                   # line number
            flush(); buf = []; continue
        if tok.endswith(']'):
            flush(); cur = dict(lemma=' '.join(buf + [tok[:-1]]), lemma_sig=None, groups=[]); buf = []; continue
        s, comma = sig_of(tok)
        if s:
            if cur is None: cur = dict(lemma=None, lemma_sig=None, groups=[])
            if not buf and cur['lemma'] is not None and not cur['groups'] and cur['lemma_sig'] is None:
                cur['lemma_sig'] = s
            else:
                cur['groups'].append((' '.join(buf), s))
            buf = []
            if not comma: flush()
            continue
        buf.append(tok)
    flush()
    return out


def units(lines, lo, hi):
    U = []
    for hl, core, partial, body in blocks(lines, lo, hi):
        txt = ' '.join(t for _, t in body)
        for e in entries(txt):
            groups = [set(g) for _, g in e['groups']]
            if e['lemma_sig']: groups = [set(e['lemma_sig'])] + groups
            allc = [w for g in groups for w in g]
            if len(allc) != len(set(allc)): continue            # witness in two states: OCR or sub-hand noise
            cited = set(allc)
            # Mommsen reports H P V L A O B X in full (prooem. p. LXXII); if OCR dropped one from the header
            # but the page cites it, it is present
            core = core | (cited & FULL)
            has_lemma = bool(e['lemma_sig'])
            if not e['lemma_sig']:
                rest = core - cited
                if rest: groups = [rest] + groups; has_lemma = True
            if len(groups) < 2: continue
            U.append(dict(page_header=hl, lemma=e['lemma'], lemma_group=0 if has_lemma else None,
                          readings=[r for r, _ in e['groups']],
                          groups=[''.join(sorted(g, key=ORDER.index)) for g in groups]))
    return U


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('txt'); ap.add_argument('--out', required=True)
    ap.add_argument('--from', dest='lo', type=int, default=1); ap.add_argument('--to', dest='hi', type=int, default=0)
    a = ap.parse_args()
    lines = open(a.txt, encoding='utf-8', errors='replace').readlines()
    hi = a.hi or len(lines)
    U = units(lines, a.lo - 1, hi)
    os.makedirs(a.out, exist_ok=True)
    json.dump(U, open(os.path.join(a.out, 'units.json'), 'w'), ensure_ascii=False, indent=0)
    wit = [w for w in ORDER if any(w in g for u in U for g in u['groups'])]
    rows = []
    for u in U:
        st = {w: str(k) for k, g in enumerate(u['groups']) for w in g}
        # informative only: at least two states held by >= 2 witnesses each
        sizes = Counter(len(g) for g in u['groups'])
        if sum(1 for g in u['groups'] if len(g) >= 2) < 2: continue
        rows.append(' '.join(st.get(w, '?') for w in wit))
    with open(os.path.join(a.out, 'jordanes.num'), 'w') as f:
        f.write(' '.join(wit) + '\n' + '\n'.join(rows) + '\n')
    pairs = Counter(u['groups'][i] for u in U for i in range(len(u['groups'])) if 2 <= len(u['groups'][i]) <= len(wit) - 2)
    print(len(U), 'units;', len(rows), 'informative; witnesses', wit)
    for g, c in pairs.most_common(30): print(f'{c:5d}  {g}')


if __name__ == '__main__':
    main()
