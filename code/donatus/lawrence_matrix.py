"""Lawrence prologue BHL 4759: build a witness x unit matrix (.num) from the apparatus of the seven copies in
mc_editions/leaves/paper/fragmentology/tools/prologue.py (APP, ROWS). Our own edition: no copyright issue.

Witnesses: T M F P B R Mi (V and S run parallel only in the last paragraphs: left out).
Rules: negative apparatus (unnamed = T's text); first hand (T^ac counts, ^pc and ^2 ignored); a segment
with sigla but no reading = those witnesses read the lemma; 'om.' = omission; spelling-only differences
(build_matrix.ortho) count as the lemma; M is '?' where its line is lost (ROWS M locus empty), where
the lemma lies in [[...]] (not legible in M), or where the apparatus says 'not in M' / 'n.l.'.
Parenthetical notes are ignored. Every unit is printed with its parse in --audit for proofing.
Usage: lawrence_matrix.py <out.num> [--audit audit.tsv]"""
import sys, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'mc_editions', 'leaves', 'paper', 'fragmentology', 'tools'))
import prologue
from build_matrix import ortho

WIT = ['T', 'M', 'F', 'P', 'B', 'R', 'Mi']
SIG = re.compile(r'^(T|M|F|P|B|R|Mi)(\^(ac|pc|2))?[,;.]?$')


def plain(s):
    s = re.sub(r'\\textsuperscript\{([^}]*)\}', r'^\1', s)
    s = re.sub(r'\\textit\{', '{', s)
    s = s.replace(r'\dots', '…').replace('\\ ', ' ')
    s = re.sub(r'\\[a-zA-Z]+', '', s)
    return s


def strip_parens(s):
    while True:
        t = re.sub(r'\([^()]*\)', '', s)
        if t == s: return s
        s = t


def lemma_text(l):
    return ' '.join(re.findall(r'\{([^{}]*)\}', plain(l))) or plain(l)


rows_by_line = {}
for loc, dip, norm, mloc in prologue.ROWS:
    rows_by_line[loc] = (norm, mloc)


def m_absent(locus, lemma):
    first = locus.split('--')[0]
    if first in rows_by_line:
        norm, mloc = rows_by_line[first]
        if not mloc: return True
        hidden = ' '.join(re.findall(r'\[\[(.*?)\]\]', norm)).lower()
        lw = [w for w in re.findall(r'[a-z]+', lemma.lower()) if len(w) > 2]
        if lw and hidden and all(w in hidden for w in lw): return True
    return False


def parse(entry):
    locus, lem, var = entry
    lem_t = lemma_text(lem)
    reads = {}       # witness -> reading text ('' = omission, None = lemma)
    absent = set()
    v = plain(var)
    if re.search(r'not in M|n\.\s*l\.', v): absent.add('M')
    for seg in re.split(r'[:;]', v):
        seg = strip_parens(seg).strip()
        if not seg: continue
        txt = ' '.join(re.findall(r'\{([^{}]*)\}', seg)).strip()
        rest = re.sub(r'\{[^{}]*\}', ' ', seg)
        toks = rest.replace(',', ' ').replace(';', ' ').split()
        sig = []
        for t in toks:
            mm = SIG.match(t)
            if mm and mm.group(3) in (None, 'ac'): sig.append(mm.group(1))
        if not sig: continue
        if 'om.' in rest: r = ''
        elif txt: r = txt
        else: r = None
        for w in sig: reads[w] = r
    if m_absent(locus, lem_t): absent.add('M')
    return locus, lem_t, reads, absent


def build():
    units = []
    for e in prologue.APP:
        locus, lem_t, reads, absent = parse(e)
        ids = {ortho(lem_t): '0'}; row = {}
        for w in WIT:
            if w in absent: row[w] = '?'; continue
            r = reads.get(w, None)
            if r is None: row[w] = '0'; continue
            key = '<om>' if r == '' else ortho(r)
            if key == ortho(lem_t): row[w] = '0'; continue
            row[w] = ids.setdefault(key, str(len(ids)))
        units.append((locus, lem_t, reads, absent, row))
    return units


if __name__ == '__main__':
    out = sys.argv[1]; aud = sys.argv[sys.argv.index('--audit') + 1] if '--audit' in sys.argv else None
    units = build()
    kept = [u for u in units if len(set(v for v in u[4].values() if v != '?')) >= 2]
    with open(out, 'w') as f:
        f.write(' '.join(WIT) + '\n')
        for u in kept: f.write(' '.join(u[4][w] for w in WIT) + '\n')
    if aud:
        with open(aud, 'w') as f:
            for locus, lem, reads, absent, row in units:
                f.write('\t'.join([locus, lem, str(reads), ','.join(sorted(absent)), ' '.join(f'{w}={row[w]}' for w in WIT)]) + '\n')
    print(f'APP entries {len(units)}  informative units {len(kept)}  M absent in {sum(u[4]["M"] == "?" for u in units)}')
