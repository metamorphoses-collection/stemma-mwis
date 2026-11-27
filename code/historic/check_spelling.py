"""CHECK (2026-11-23): orthographic classification of extracted readings for the series test.
Not part of the frozen pipeline. A departure (one group's reading in one unit) is 'spelling-only' when, after an
orthographic normalisation, it equals its lemma but differs from it as written.

Lemma: the unit's own lemma where the apparatus prints one ('word]'); otherwise recovered from the edited text of the
same leaf: the reading is spelling-only if some word of that leaf's edited text (lines the frozen extractor does not
read as apparatus) normalises to the same form but is written differently, and the reading itself does not occur
verbatim in that text. Multi-word readings must match word for word (consecutive words of the text).
Readings with editorial words (om., add., del., sup. lin., ...), digits or sigla are never spelling-only.

Two normalisations:
  strict: u/v, i/j/y, k/c, ae/oe -> e, h dropped (ch, ph->f, th), ti/ci before a vowel, mpn/mpt/mps, assimilation of
          ad-/in-/con-/sub-/ob- prefixes (np/nb/nm, dl/dn/dp/dr/ds/dt/dc/df/dg, nl/nr, bp/bm/bc), bt/bs -> pt/ps, xs -> x,
          quu -> cu, final d -> t, doubled letters single;
  broad:  strict + e/i and o/u merged, b/u (betacism), final m dropped  (Late Latin phonetic spelling; an upper bound
          that also merges some grammatical endings, e.g. natale/natali, videt/vidit).
"""
import re, os, json, sys, unicodedata
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mgh_rooted_extract import make, is_app

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'bedier')
EDIT = re.compile(r"^(om|add|del|eras|corr|sup|lin|mg|marg|ras|post|ante|cf|coni|man|in|ex|s|m|al|uel|vel|scr|ss|exp|praem|transp|sec|pr|alt|ead|eadem|rec|tit|inscr|subscr|lac|lacuna|deest|desunt|des|fol|cod|codd|edd|ed|recte|sic|bis|iter|alterum|alt)\.?$")


def _deacc(s): return ''.join(c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn')


def norm(w, broad=False):
    w = _deacc(w.lower()).replace('*', '')
    w = re.sub(r'[^a-z]', '', w)
    w = w.replace('v', 'u').replace('j', 'i').replace('y', 'i').replace('k', 'c')
    w = w.replace('ae', 'e').replace('oe', 'e').replace('ph', 'f')
    w = w.replace('h', '')
    w = re.sub(r'ci(?=[aeiou])', 'ti', w)
    w = w.replace('mpn', 'mn').replace('mpt', 'mt').replace('mps', 'ms')
    w = re.sub(r'^(a)d(?=[clnprstfg])', lambda m: 'a' + '#', w)            # ad- assimilation: adp- ~ app-
    w = re.sub(r'a#(.)', r'a\1\1', w)
    w = re.sub(r'^(i|co)n(?=[lr])', lambda m: m.group(1) + '#', w); w = re.sub(r'(i|co)#(.)', r'\1\2\2', w)
    w = re.sub(r'n(?=[pbm])', 'm', w)
    w = re.sub(r'^(su|o)b(?=[pmcfg])', lambda m: m.group(1) + '#', w); w = re.sub(r'(su|o)#(.)', r'\1\2\2', w)
    w = w.replace('bt', 'pt').replace('bs', 'ps').replace('xs', 'x').replace('quu', 'cu')
    w = re.sub(r'd$', 't', w)
    if broad:
        w = w.replace('e', 'i').replace('o', 'u').replace('b', 'u')
        w = re.sub(r'm$', '', w)
    return re.sub(r'(.)\1+', r'\1', w)


def words(s): return [t for t in re.split(r'\s+', s.strip()) if t]


def clean_reading(r):
    r = re.sub(r'\([^)]*\)', ' ', r or '')
    return [t.strip(',.;:!?|') for t in words(r)]


ROMAN = re.compile(r'^m{0,4}(cm|cd|d?c{0,4})(xc|xl|l?x{0,4})(ix|iu|iv|[uv]?i{0,4})$')
NOT_NUM = {'mi', 'di', 'ui', 'vi', 'dic', 'lux', 'mix', 'uix', 'vix', 'mil', 'cui', 'xi'}


def is_numeral(t):
    t = t.lower()
    return len(t) >= 2 and t not in NOT_NUM and bool(ROMAN.match(t))


def eligible(r, sigla=''):
    toks = clean_reading(r)
    if not toks: return False
    for t in toks:
        if EDIT.match(t.lower()) or re.search(r'\d', t) or not re.fullmatch(r"[A-Za-zÀ-ÿ*'\-]+", t): return False
        if t.isupper() and len(t) <= 6: return False                                     # sigla or editorial capitals
        if '*' in t or is_numeral(t): return False                                        # erasures, numerals
    return True


def spelling_vs(reading, lemma, broad):
    a, b = clean_reading(reading), clean_reading(lemma)
    if len(a) != len(b) or not a: return False
    na, nb = [norm(x, broad) for x in a], [norm(x, broad) for x in b]
    return na == nb and [x.lower() for x in a] != [x.lower() for x in b] and all(na)


class Leaves:
    """Edited text of a text's leaves, from the same vision transcription the units were extracted from, with printed
    line numbers estimated from the marginal numbers (multiples of 5, at either end of a line). A unit's lemma window is
    the text lines numbered N-1..N+1, N being the line number that opens its apparatus entry; if N or the numbering
    cannot be found, the whole leaf is the window."""
    def __init__(self, text_id, settings):
        lv = json.load(open(os.path.join(ROOT, 'runs', text_id, 'leaves.json')))
        sig_of = make(settings['sigla'], settings['ocr'], settings['digits'], settings.get('lower', ''))
        buf = ''
        for l in range(lv['leaves'][0], lv['leaves'][1] + 1):
            f = os.path.join(ROOT, 'ocr_llm', lv['vol'], f"{lv['vol']}_{l:05d}.llm.txt")
            if os.path.exists(f): buf += f'@@LEAF {l}\n' + open(f, encoding='utf-8').read() + '\n'
        self.L = L = buf.split('\n')
        self.text = {}                        # leaf -> list of (est_line_no or None, line text)
        self.found_N = self.window_used = 0
        leaf = None; start = None; rows = []
        def close_leaf():
            if leaf is None: return
            anchors = []
            for k, (_, t) in enumerate(rows):
                m = re.match(r'^(\d{1,3})\s+\S', t) or re.search(r'\S\s+(\d{1,3})$', t)
                if m and int(m.group(1)) % 5 == 0: anchors.append((k, int(m.group(1))))
            est = []
            for k, (_, t) in enumerate(rows):
                c = [n + (k - a) for a, n in anchors if abs(k - a) <= 4]
                est.append(set(c))
            self.text[leaf] = [(est[k], t) for k, (_, t) in enumerate(rows)]
        for i, raw in enumerate(L):
            t = raw.strip()
            m = re.match(r'^@@LEAF (\d+)', t)
            if m: close_leaf(); leaf = int(m.group(1)); rows = []; start = None; continue
            if is_app(t, sig_of, in_block=start is not None):
                if start is None: start = i + 1
                continue
            if t:
                nxt = [L[j].strip() for j in range(i + 1, min(i + 6, len(L))) if L[j].strip()][:2]
                if start is not None and any(is_app(x, sig_of, True) for x in nxt): continue
                start = None
                rows.append((i + 1, t))
        close_leaf()

    def entry_number(self, u, reading):
        """printed line number that opens the unit's apparatus entry (None if not found)."""
        el = u.get('entry_line') or u.get('line')
        if not el: return None
        first = (clean_reading(reading) or [''])[0]
        for ln in (el, el + 1):
            if ln - 1 >= len(self.L): break
            t = self.L[ln - 1]
            p = t.find(first) if first else -1
            if p < 0: continue
            nums = re.findall(r'(?:^|\s)(\d{1,3})(?=[\s.,]|$)', t[:p])
            if nums: return int(nums[-1])
            prev = re.findall(r'(?:^|\s)(\d{1,3})(?=[\s.,]|$)', self.L[ln - 2]) if ln >= 2 else []
            return int(prev[-1]) if prev else None
        return None

    def window(self, u, reading):
        rows = self.text.get(u['leaf']) or []
        N = self.entry_number(u, reading)
        if N is not None:
            sel = [t for est, t in rows if est & {N - 1, N, N + 1}]
            if sel:
                self.window_used += 1
                return sel
        return [t for _, t in rows]

    def spelling(self, u, reading, broad):
        txt = re.sub(r'(\w)-\s+(\w)', r'\1\2', ' '.join(self.window(u, reading)))
        W = [t.strip(',.;:!?|()[]"\'') for t in words(txt)]
        W = [w for w in W if w]
        a = clean_reading(reading)
        if not a or not W: return False
        k = len(a); na = [norm(x, broad) for x in a]; la = [x.lower() for x in a]
        hit_norm = hit_exact = False
        for i in range(len(W) - k + 1):
            seg = W[i:i + k]; ls = [x.lower() for x in seg]
            if ls == la: hit_exact = True; break
            if [norm(x, broad) for x in seg] == na: hit_norm = True
        return hit_norm and not hit_exact and all(na)


def classify_unit(u, leaves, broad, union=False):
    """list of booleans, one per departure group (groups[1:]): True = spelling-only. Lemma printed in the apparatus ->
    compared with it; no lemma -> the line window of the edited text. union=True: spelling-only if either says so."""
    out = []
    for r in u['readings']:
        if not eligible(r): out.append(False); continue
        lem = spelling_vs(r, u['lemma'], broad) if u.get('lemma') and eligible(u['lemma']) else None
        if union: out.append(bool(lem) or leaves.spelling(u, r, broad))
        elif lem is not None: out.append(lem)
        elif u.get('lemma'): out.append(False)
        else: out.append(leaves.spelling(u, r, broad))
    return out
