"""A6 pipeline (PREREG_BEDIER_20261113.md): per-leaf multi-engine OCR -> voted units.
For each text: printed page range (first..last page in its IA range, as in bedier_settings.pages_in) -> leaves whose printed
page number (read on the leaf's own OCR head) falls in that range; each variant (strip2, strip3, tess, tess4) is extracted
with the text's frozen settings over those leaves (@@LEAF markers); on each leaf a group is kept with the multiplicity at
which at least two of the four variants report it (any group size). Each kept occurrence becomes a one-group unit.
Usage: bedier_a6.py OCR_ROOT --out RUNS_A6 [--only text_id,...]
"""
import json, os, re, sys, glob, subprocess, argparse, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
B = os.path.expanduser('~/MC/projects/priscus_excerpta/data/historic/bedier')
VARIANTS = ['strip2', 'strip3', 'tess', 'tess4']
VOLS = {'AA 1': ['mgh_auct_ant_1_1', 'mgh_auct_ant_1_2'], 'AA 2': ['mgh_auct_ant_2'], 'AA 3': ['mgh_auct_ant_3_1', 'mgh_auct_ant_3_2'],
        'AA 4': ['mgh_auct_ant_4_1', 'mgh_auct_ant_4_2'], 'AA 5': ['mgh_auct_ant_5_1', 'mgh_auct_ant_5_2'], 'AA 6.1': ['mgh_auct_ant_6_1'],
        'AA 6.2': ['mgh_auct_ant_6_2'], 'AA 7': ['mgh_auct_ant_7'], 'AA 8': ['mgh_auct_ant_8'], 'AA 9': ['mgh_auct_ant_9'],
        'AA 10': ['mgh_auct_ant_10'], 'AA 11': ['mgh_auct_ant_11'], 'AA 12': ['mgh_auct_ant_12'], 'AA 13': ['mgh_auct_ant_13'],
        'AA 14': ['mgh_auct_ant_14'], 'AA 15': ['mgh_auct_ant_15']}
OVERRIDE_START = {'prosper_epitoma_chronicon': ('mgh_auct_ant_9', 357, 0.75)}   # A8c: start match 0.749 at leaf 357 (p. 346, Prosper preface), admitted by CJ
MIN_LEAVES = 50          # a volume counts as OCRed only when its phase-2 output is present (phase 1 left 3 leaves in some)


def page_range(s):
    import importlib.util
    spec = importlib.util.spec_from_file_location('bs', os.path.join(HERE, 'bedier_settings.py'))
    src = open(os.path.join(HERE, 'bedier_settings.py')).read().split('\nout = []')[0]
    ns = {}; exec(src, ns)
    lines = open(os.path.join(B, 'ocr', s['ia_id'] + '.txt'), encoding='utf-8', errors='replace').readlines()
    P = ns['pages_in'](lines, *s['lines'])
    return (P[0], P[-1]) if P else None


def leaf_pages(voldir, vol):
    """printed page number of each leaf, read from the first lines of its strip2 (else tess) OCR."""
    out = {}
    for f in sorted(glob.glob(os.path.join(voldir, f'{vol}_*.strip2.txt'))):
        leaf = int(re.search(r'_(\d{5})\.strip2', f).group(1))
        head = [l.strip() for l in open(f, encoding='utf-8', errors='replace').read().split('\n')[:4] if l.strip()]
        for t in head:
            m = re.fullmatch(r'(\d{1,3})', t) or re.match(r'^(\d{1,3})\s+\D{4,}$', t) or re.match(r'^\D{4,}\s+(\d{1,3})$', t)
            if m: out[leaf] = int(m.group(1)); break
    # keep only readings consistent with a locally constant offset (majority offset among neighbours)
    off = {l: l - p for l, p in out.items()}
    good = {}
    L = sorted(off)
    for i, l in enumerate(L):
        nb = [off[x] for x in L[max(0, i - 3):i + 4]]
        if collections.Counter(nb).most_common(1)[0][0] == off[l]: good[l] = out[l]
    return good


import difflib
norm = lambda x: re.sub(r'[^a-z]', '', x.lower())


def anchor_lines(lines, i, step, k=4, maxscan=600):
    """k lines of edited text (>= 6 words, >= 80% lower-case alphabetic words, no sigla runs) from IA line i in direction step."""
    out = []; j = i; n = 0
    while 0 <= j < len(lines) and len(out) < k and n < maxscan:
        t = lines[j].strip(); w = t.split()
        if len(w) >= 6 and sum(1 for x in w if re.fullmatch(r"[a-z][a-z,.;:!?']*", x)) >= 0.8 * len(w) and not re.search(r'\b[A-Z]{2,}\b', t):
            out.append(t)
        j += step; n += 1
    return out


def locate_leaf(roots, anchors, after=None):
    """(vol, leaf, ratio) of the leaf whose strip2 OCR best matches the anchor lines."""
    best = (0.0, None, None); NA = [norm(a) for a in anchors]
    for vd, vol in roots:
        for f in glob.glob(os.path.join(vd, f'{vol}_*.strip2.txt')):
            lf = int(re.search(r'_(\d{5})\.strip2', f).group(1))
            if after and (vol != after[0] or lf < after[1]): continue
            L = [norm(x) for x in open(f, encoding='utf-8', errors='replace').read().split('\n')]
            L = [x for x in L if len(x) >= 15]
            sc = 0.0
            for na in NA:
                cand = [x for x in L if abs(len(x) - len(na)) < 0.5 * len(na)]
                sc += max((difflib.SequenceMatcher(None, na, x).ratio() for x in cand), default=0.0)
            sc /= max(1, len(NA))
            if sc > best[0]: best = (sc, vol, int(re.search(r'_(\d{5})\.strip2', f).group(1)))
    return best


def head_of(vd, vol, leaf):
    f = os.path.join(vd, f'{vol}_{leaf:05d}.strip2.txt')
    if not os.path.exists(f): return ''
    for t in open(f, encoding='utf-8', errors='replace').read().split('\n')[:3]:
        n = norm(t)
        if len(n) >= 8: return n
    return ''


def extend_end(vd, vol, l1, cap=15):
    """A8: after the content-located end leaf, add the following leaves whose running head (first line of >= 8 letters,
    letters only) matches the end leaf's head (difflib ratio >= 0.8); stop at the first leaf that does not, at most cap."""
    h = head_of(vd, vol, l1); l = l1
    if not h: return l1
    while l - l1 < cap:
        hn = head_of(vd, vol, l + 1)
        if not hn or difflib.SequenceMatcher(None, h, hn).ratio() < 0.8: break
        l += 1
    return l


def toc_leaves(vol):
    """A8b: leaves where a dMGH table-of-contents section starts (page labels only; titles are not used)."""
    import urllib.request
    d = os.path.join(B, 'dmgh_toc'); os.makedirs(d, exist_ok=True); f = os.path.join(d, vol + '.json')
    if not os.path.exists(f):
        js = urllib.request.urlopen(f'https://www.dmgh.de/{vol}/book.js', timeout=60).read().decode('utf-8', 'replace')
        pages = json.loads(re.search(r'var pages = (\[.*?\]);', js, re.S).group(1))
        toc = js[js.index('br.updateTOC('):]
        labels = re.findall(r'"pagenum":\s*"([^"]*)"', toc)
        json.dump(dict(labels=labels, leaves=sorted(pages.index(p) + 1 for p in labels if p in pages), n_leaves=len(pages)), open(f, 'w'))
    return json.load(open(f))['leaves']


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('root'); ap.add_argument('--out', required=True); ap.add_argument('--only', default='')
    ap.add_argument('--source', choices=['vote', 'llm'], default='vote'); ap.add_argument('--llm-root', default=os.path.join(B, 'ocr_llm'))
    ap.add_argument('--list-leaves', action='store_true')
    a = ap.parse_args()
    st = json.load(open(os.path.join(B, 'settings.json')))
    only = set(filter(None, a.only.split(',')))
    for s in st:
        t = s['text_id']
        if only and t not in only: continue
        if len(s['sigla']) < 2: continue
        roots = [(os.path.join(a.root, v), v) for v in VOLS[s['volume']]]
        roots = [(d, v) for d, v in roots if os.path.isdir(d) and len(glob.glob(os.path.join(d, '*.strip2.txt'))) >= MIN_LEAVES]
        if not roots: print(t, 'volume not yet OCRed:', VOLS[s['volume']]); continue
        ia = open(os.path.join(B, 'ocr', s['ia_id'] + '.txt'), encoding='utf-8', errors='replace').readlines()
        span = s['lines'][1] - s['lines'][0]
        a0 = anchor_lines(ia, s['lines'][0] + 5, +1) or anchor_lines(ia, s['lines'][0] + 5, +1, maxscan=span)   # A8: widen only if none
        a1 = anchor_lines(ia, s['lines'][1] - 1, -1) or anchor_lines(ia, s['lines'][1] - 1, -1, maxscan=span)
        r0, v0, l0 = locate_leaf(roots, a0)
        r1, v1, l1 = locate_leaf(roots, a1, after=(v0, l0))
        end_by = 'content'
        if t in OVERRIDE_START:                                         # A8c: CJ's decision, disclosed as a deviation
            v0, l0, r0 = OVERRIDE_START[t]
        if v0 is not None and r0 >= 0.75 and (v1 != v0 or r1 < 0.75 or l1 is None or l1 < l0):
            nxt = [x for x in toc_leaves(v0) if x > l0]                 # A8b: end at the next dMGH section boundary
            if nxt: v1, l1, r1, end_by = v0, nxt[0] - 1, 1.0, 'toc'
        if v0 is None or v0 != v1 or r0 < 0.75 or r1 < 0.75 or l1 < l0:
            print(t, f'not located: start {v0} {l0} ({r0:.2f}) end {v1} {l1} ({r1:.2f})'); continue
        vol = v0; vd = os.path.join(a.root, vol)
        l1x = extend_end(vd, vol, l1) if end_by == 'content' else l1                                   # A8: continuation leaves with the same running head
        leaves = list(range(l0, l1x + 1)); pr = (l0, l1x)
        od = os.path.join(a.out, t); os.makedirs(od, exist_ok=True)
        if a.source == 'llm':
            txt = os.path.join(od, 'llm.txt'); missing = []
            with open(txt, 'w') as fh:
                for l in leaves:
                    f = os.path.join(a.llm_root, vol, f'{vol}_{l:05d}.llm.txt')
                    if os.path.exists(f): fh.write(f'@@LEAF {l}\n' + open(f, encoding='utf-8').read() + '\n')
                    else: missing.append(l)
            if missing: print(t, 'LLM reading missing for leaves', missing[:10], len(missing)); continue
            cmd = ['python3', os.path.join(HERE, 'mgh_rooted_extract.py'), txt, '--sigla', s['sigla'], '--out', os.path.join(od, 'llm'), '--name', t, '--digits', s['digits'], '--lower', s.get('lower', '')]
            if s['ocr']: cmd += ['--ocr', s['ocr']]
            subprocess.run(cmd, capture_output=True, check=True)
            units = json.load(open(os.path.join(od, 'llm', 'units.json')))
            json.dump(units, open(os.path.join(od, 'units.json'), 'w'))
            json.dump(dict(vol=vol, leaves=[l0, l1x], located_end=l1, end_by=end_by, match=[round(r0, 2), round(r1, 2)], source='llm'), open(os.path.join(od, 'leaves.json'), 'w'))
            print(f"{t:32s} {vol:18s} leaves {l0}-{l1x} ({len(leaves)}) units {len(units)} (shared {sum(1 for u in units if any(len(g) >= 2 for g in u['groups'][1:]))})")
            continue
        if a.list_leaves:
            print('LEAVES', t, vol, l0, l1x, end_by); continue
        per = {}
        for v in VARIANTS:
            txt = os.path.join(od, f'{v}.txt')
            with open(txt, 'w') as fh:
                for l in leaves:
                    f = os.path.join(vd, f'{vol}_{l:05d}.{v}.txt')
                    if os.path.exists(f): fh.write(f'@@LEAF {l}\n' + open(f, encoding='utf-8', errors='replace').read() + '\n')
            cmd = ['python3', os.path.join(HERE, 'mgh_rooted_extract.py'), txt, '--sigla', s['sigla'], '--out', os.path.join(od, v), '--name', t, '--digits', s['digits']]
            if s['ocr']: cmd += ['--ocr', s['ocr']]
            subprocess.run(cmd, capture_output=True, check=True)
            U = json.load(open(os.path.join(od, v, 'units.json')))
            c = collections.defaultdict(collections.Counter)
            for u in U:
                for g in u['groups'][1:]:
                    if g: c[u['leaf']][''.join(sorted(g))] += 1
            per[v] = c
        units = []
        for l in leaves:
            keys = set().union(*[per[v][l] for v in VARIANTS])
            for k in keys:
                m = sorted((per[v][l][k] for v in VARIANTS), reverse=True)[1]
                g = ''.join(x for x in s['sigla'] if x in k)          # witness order of the settings
                units += [dict(line=None, entry_line=None, leaf=l, lemma=None, lemma_group=0, groups=['', g], readings=[''])] * m
        json.dump(units, open(os.path.join(od, 'units.json'), 'w'))
        json.dump(dict(vol=vol, leaves=[l0, l1x], located_end=l1, match=[round(r0, 2), round(r1, 2)]), open(os.path.join(od, 'leaves.json'), 'w'))
        print(f"{t:32s} {vol:18s} leaves {l0}-{l1} ({len(leaves)}) match {r0:.2f}/{r1:.2f} units {len(units)} (shared {sum(1 for u in units if len(u['groups'][1]) >= 2)})")


if __name__ == '__main__':
    main()
