"""CHECK (2026-11-23): Weitzman comparison by editor. Applies the one counting rule of
the counting of two-way tops (supplementary notes) to the blind top records of MGH AA 1-15 and SRM 1-7: explicit if
filiation stated in a diagram or in words, confidence high or medium, not flagged classes_only / recensions /
editor_undecided / several_archetypes, and k_E >= 2. Reports explicit tops per editor, the two-way share per editor,
and editor-level (one vote per editor) statistics against Weitzman's 0.71 (Latin) and 0.77 (Greek).
Records: AA = data/bedier/tops (this repo); SRM = data/bedier_srm/tops (this repo).
Two sets: 'paper' (default; the paper's 94 = AA tops with classification_A3 == explicit, the verified blind step, + SRM
by the rule) and 'explicit' (the rule applied literally to the raw AA agent flags as well: 34 AA).
Usage: check_weitzman.py [paper|explicit]
Output: data/bedier/checks_20261123/weitzman_by_editor_<set>.json"""
import json, glob, os, re, collections, ast
from math import comb
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'bedier')
SRM = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'bedier_srm', 'tops')
BAD = {'classes_only', 'recensions', 'editor_undecided', 'several_archetypes'}


def lst(x):
    if isinstance(x, list): return x
    try: return ast.literal_eval(x)
    except Exception: return []


def editor_key(e):
    """surname, for grouping 'Th. Mommsen', 'Theodor Mommsen', 'Theodorus Mommsen' etc."""
    e = re.sub(r'\(.*?\)', '', e or '?')
    e = e.split(',')[0].split(' and ')[0].split(';')[0].split('/')[0].strip()
    return e.split()[-1] if e.split() else '?'


rows = []
for series, pat in (('AA', os.path.join(R, 'tops', '*.json')), ('SRM', os.path.join(SRM, '*.json'))):
    for f in sorted(glob.glob(pat)):
        for t in json.load(open(f))['texts']:
            fl = set(lst(t.get('flags', [])))
            k = len(lst(t.get('top_branches', [])))
            ok = t.get('filiation_stated') in ('diagram', 'verbal') and t.get('confidence') in ('high', 'medium') and not (fl & BAD) and k >= 2
            paper = (t.get('classification_A3') == 'explicit') if series == 'AA' else ok     # the paper's 94: AA by A3, SRM by the rule
            rows.append(dict(series=series, text=t['text_id'], editor=t.get('editor'), ed=editor_key(t.get('editor')), k=k, explicit=ok, paper=paper))
import sys
SET = sys.argv[1] if len(sys.argv) > 1 else 'paper'
ex = [r for r in rows if r[SET]]
by = collections.defaultdict(list)
for r in ex: by[r['ed']].append(r['k'] == 2)
per = sorted(((e, len(v), sum(v)) for e, v in by.items()), key=lambda x: -x[1])
n_ed = len(per); share_mean = sum(s / n for _, n, s in per) / n_ed
all2 = sum(1 for _, n, s in per if s == n); maj2 = sum(1 for _, n, s in per if 2 * s > n)
def binom_ge(k, n, p): return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))
out = dict(set=SET, n_recorded=len(rows), n_explicit=len(ex), by_series={s: dict(explicit=sum(1 for r in ex if r['series'] == s), two_way=sum(1 for r in ex if r['series'] == s and r['k'] == 2)) for s in ('AA', 'SRM')},
           two_way=sum(1 for r in ex if r['k'] == 2), per_editor=[dict(editor=e, explicit=n, two_way=s, names=sorted({r['editor'] for r in ex if r['ed'] == e})) for e, n, s in per],
           n_editors=n_ed, editor_mean_two_way_share=share_mean, editors_all_two_way=all2, editors_majority_two_way=maj2,
           p_editor_vote_vs_weitzman={str(p): binom_ge(maj2, n_ed, p) for p in (0.71, 0.77)},
           p_text_level_vs_weitzman={str(p): binom_ge(sum(1 for r in ex if r['k'] == 2), len(ex), p) for p in (0.71, 0.77)},
           multi_branch=[r for r in ex if r['k'] != 2])
json.dump(out, open(os.path.join(R, 'checks_20261123', f'weitzman_by_editor_{SET}.json'), 'w'), indent=1, ensure_ascii=False)
print(json.dumps({k: v for k, v in out.items() if k not in ('per_editor',)}, indent=1, ensure_ascii=False))
for d in out['per_editor']: print(d)
