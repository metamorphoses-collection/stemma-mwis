"""Passio Donati III (BHL 2294): reader transcriptions (JSON per page) -> MWIS matrix (.num).

The transcriptions and the matrix are derived from Licciardello 2018 (in copyright) and stay in
a local, non-public folder; this script holds no apparatus content, only the stemma's group definitions
(p. 172) and the conversion rules.

Rules
- Columns: the 54 sigla of the stemma (p. 172); editions (BUON, DRAG, LAZZ, MOMB) ignored.
- Negative apparatus: a witness not named reads the text (0) unless absent there ([om. X] notes).
- Group sigla expand to members (closure: nu includes pi..upsilon, mu includes xi, omicron, phi;
  xi includes phi), minus exceptions. --no-expand: members named only through a group siglum get '?'.
- Hands: first hand (^ac) by default; --hand pc uses the corrected state.
- Spelling-only readings count as the lemma; units with no non-spelling reading are dropped.
- A witness assigned to two readings in one unit -> '?'.
Usage: build_matrix.py <transcr_dir> <out.num> [--no-expand] [--hand ac|pc] [--manifest out.tsv]
"""
import sys, json, glob, re, os, argparse

O2 = ['Gb', 'L1', 'N2', 'Nn', 'P1', 'V2']
BASE = {
    'kappa': ['H', 'Li', 'Me', 'Z'],
    'lambda': ['B2', 'Ko'],
    'phi': ['C1', 'C2', 'D', 'P2', 'SO'],
    'xi': ['P', 'U1', 'U2', 'U3', '@phi'],
    'omicron': ['Bo', 'F16', 'Ri'],
    'mu': ['Ar', '@xi', '@omicron'],
    'pi': ['F17', 'Pg'],
    'rho': ['F15', 'Fa', 'T1'],
    'sigma': ['Ad', 'L3', 'R4'],
    'tau': ['F14', 'F19', 'L2', 'L4', 'Pi', 'R3', 'V3'],
    'upsilon': ['F20', 'L5', 'Mi', 'T2'],
    'nu': ['Au', 'F18', 'L6', 'M1', 'M2', 'N1', 'P3', 'P4', 'R2', 'Re',
           '@pi', '@rho', '@sigma', '@tau', '@upsilon'],
}
EDITIONS = {'BUON', 'DRAG', 'LAZZ', 'MOMB'}


def members(g):
    out = []
    for m in BASE[g]:
        out += members(m[1:]) if m.startswith('@') else [m]
    return out


GROUPS = {g: members(g) for g in BASE}
WIT = O2 + [w for g in ['kappa', 'lambda', 'mu', 'nu'] for w in GROUPS[g]]
assert len(WIT) == len(set(WIT)) == 54, len(WIT)
ALIAS = {'tutti': 'ALL', 'κ': 'kappa', 'x': 'kappa', 'λ': 'lambda', 'μ': 'mu', 'ν': 'nu', 'v': 'nu', 'ξ': 'xi',
         'ο': 'omicron', 'o': 'omicron', 'π': 'pi', 'ρ': 'rho', 'ϱ': 'rho', 'σ': 'sigma', 'τ': 'tau',
         'υ': 'upsilon', 'φ': 'phi', 'ϕ': 'phi'}


def ortho(t):
    """Orthographic skeleton: spelling carries no stemmatic weight (Licciardello p. 162)."""
    t = t.lower().replace('ę', 'e').replace('æ', 'e').replace('ae', 'e').replace('oe', 'e')
    t = t.replace('j', 'i').replace('y', 'i').replace('v', 'u').replace('k', 'c')
    t = re.sub(r'ti(?=[aeiou])', 'ci', t)
    t = t.replace('michi', 'mihi').replace('nichil', 'nihil')
    t = t.replace('mpn', 'mn').replace('ph', 'f').replace('th', 't').replace('ch', 'c').replace('h', '')
    t = re.sub(r'(?<=[aeiou])m(?=[bp])', 'n', t)          # inm-/imm-, conp-/comp- etc.
    t = re.sub(r'\bin(?=[mbp])', 'im', t); t = re.sub(r'\bcon(?=[mbp])', 'com', t)
    t = re.sub(r'\bad(?=[ftcpgls])', 'a', t)              # adf-/aff- assimilation
    t = re.sub(r'([a-z])\1+', r'\1', t)                   # double consonants
    t = re.sub(r'[^a-z ]', '', t)
    return ' '.join(t.split())


def norm(s):
    s = s.strip()
    return ALIAS.get(s, s)


def parse_token(tok):
    """-> (set of witnesses, via_group: bool, hand or None, unknown token or None)"""
    tok = tok.strip()
    hand = None
    m = re.match(r'^(.*)\^(ac|pc)$', tok)
    if m:
        tok, hand = m.group(1), m.group(2)
    m = re.match(r'^([^(]+)\((.*)\)$', tok)
    exc = []
    if m:
        tok = m.group(1)
        exc = [norm(re.sub(r'^(-|om[:.]?\s*)', '', e.strip())) for e in m.group(2).split(',') if e.strip()]
    tok = norm(tok)
    if tok in EDITIONS:
        return set(), False, hand, None
    if tok in GROUPS or tok == 'ALL':
        s = set(WIT) if tok == 'ALL' else set(GROUPS[tok])
        for e in exc:
            s -= set(GROUPS[e]) if e in GROUPS else {e}
        return s, True, hand, None
    if tok in WIT:
        return {tok}, False, hand, None
    return set(), False, hand, tok


def rejoin(toks):
    """Readers sometimes split 'nu(-F19,-N1)' at the commas; glue unbalanced pieces back."""
    out, buf = [], ''
    for t in toks:
        buf = f'{buf},{t}' if buf else t
        if buf.count('(') <= buf.count(')'):
            out.append(buf.replace(' ', '')); buf = ''
    if buf: out.append(buf)
    return out


def om_exceptions(tok):
    """'sigma(om:Ad)' = the group without Ad, and Ad LACKS the passage (not 'tranne')."""
    m = re.match(r'^[^(]+\((.*)\)', tok.replace(' ', ''))
    out = set()
    if m:
        for e in m.group(1).split(','):
            mm = re.match(r'^om[:.]?(.+)$', e)
            if mm:
                x = norm(mm.group(1)); out |= set(GROUPS[x]) if x in GROUPS else {x}
    return out


def load(tdir):
    pages = sorted(glob.glob(os.path.join(tdir, 'p*.json')))
    return [json.load(open(p)) for p in pages]


def build(pages, expand=True, hand='ac'):
    absent = {}
    for pg in pages:
        for c in pg.get('coverage', []):
            k = (str(c.get('chapter')), str(c.get('sentence')))
            for t in rejoin(c.get('wits', [])):
                s, _, _, _ = parse_token(t)
                absent.setdefault(k, set()).update(s)
    rows, manifest, unknown = [], [], {}
    stats = dict(units=0, spelling_only=0, kept=0, conflicts=0, group_tokens=0)
    for pg in pages:
        for u in pg.get('units', []):
            stats['units'] += 1
            k = (str(u.get('chapter')), str(u.get('sentence')))
            lo = ortho(u.get('lemma', ''))
            reads = [r for r in u.get('readings', []) if r.get('type') != 'spelling'
                     and not (ortho(r.get('text', '')) == lo and r.get('type') not in ('omission', 'order'))]
            if not reads:
                stats['spelling_only'] += 1
                continue
            cell = {w: '0' for w in WIT}
            for w in absent.get(k, ()):
                cell[w] = '?'
            assigned = {}
            unsure = set()
            other_hand = set()
            gone = set()
            lemma_like = []
            for r in u.get('readings', []):
                if r.get('type') == 'spelling':
                    lemma_like.append(r)
            ids = {ortho(u.get('lemma', '')): '0'}
            for r in u.get('readings', []):
                t = ortho(r.get('text', ''))
                if r.get('type') == 'spelling' or (ids.get(t) == '0' and r.get('type') not in ('omission', 'order')):
                    v = '0'
                elif r.get('type') == 'omission':
                    v = ids.setdefault('<om>', str(len(ids)))
                else:
                    v = ids.setdefault(t, str(len(ids)))
                nonreading = r.get('text', '').strip().startswith('[')  # [illeggibile], [eraso]: state unknown
                for t in rejoin(r.get('wits', [])):
                    s, via, h, unk = parse_token(t)
                    if nonreading:
                        if not h or h == hand:
                            gone |= s
                        continue
                    gone |= om_exceptions(t)
                    if unk:
                        unknown[unk] = unknown.get(unk, 0) + 1
                    if h and h != hand:
                        other_hand |= s   # this witness's state under the chosen hand is unknown unless listed
                        continue
                    if via:
                        stats['group_tokens'] += 1
                        if not expand:
                            unsure |= s
                            continue
                    for w in s:
                        assigned.setdefault(w, set()).add(v)
            for w in unsure | other_hand:
                if w not in assigned:
                    cell[w] = '?'
            for w, vs in assigned.items():
                if cell[w] == '?' and k in absent and w in absent[k]:
                    continue  # absent there; an explicit listing would contradict the coverage note
                if len(vs) > 1:
                    cell[w] = '?'
                    stats['conflicts'] += 1
                else:
                    cell[w] = next(iter(vs))
            for w in gone:
                cell[w] = '?'
            row = [cell[w] for w in WIT]
            if len(set(x for x in row if x != '?')) < 2:
                continue
            rows.append(row)
            manifest.append((pg['page'], u.get('chapter'), u.get('sentence'), len(reads)))
            stats['kept'] += 1
    return rows, manifest, stats, unknown, absent


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('tdir'); ap.add_argument('out')
    ap.add_argument('--no-expand', action='store_true'); ap.add_argument('--hand', default='ac')
    ap.add_argument('--manifest')
    a = ap.parse_args()
    rows, man, st, unk, absent = build(load(a.tdir), expand=not a.no_expand, hand=a.hand)
    with open(a.out, 'w') as f:
        f.write(' '.join(WIT) + '\n')
        for r in rows:
            f.write(' '.join(r) + '\n')
    if a.manifest:
        with open(a.manifest, 'w') as f:
            for m in man:
                f.write('\t'.join(map(str, m)) + '\n')
    miss = {w: sum(r[i] == '?' for r in rows) for i, w in enumerate(WIT)}
    print(json.dumps(dict(stats=st, unknown_tokens=unk, rows=len(rows),
                          missing_per_witness={w: c for w, c in miss.items() if c})))
