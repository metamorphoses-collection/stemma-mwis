"""CHECK (2026-11-23): validation of check_spelling.py. (i) On units whose apparatus prints a lemma, agreement between
the lemma-based decision and the line-window decision (which is what lemma-less units rely on); (ii) a random sample
of shared departures (>= 2 analysis witnesses) classified spelling-only under each normalisation, with the matched
text word, for reading by eye. Output: data/bedier/checks_20261123/spelling_check.json"""
import json, os, sys, random, collections, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_spelling as RS
from bedier_sensitivity import TEXTS, ROOT
S = {s['text_id']: s for s in json.load(open(os.path.join(ROOT, 'settings.json')))}
out = dict(agreement={}, sample={}, lemma_coverage={})
samp = collections.defaultdict(list)
for t in TEXTS:
    d = json.load(open(os.path.join(ROOT, 'classification_a10', t + '.json'))); Ws = set(d['analysis_witnesses'])
    U = json.load(open(os.path.join(ROOT, 'runs', t, 'units.json'))); L = RS.Leaves(t, S[t])
    ag = {}
    for broad in (False, True):
        cm = collections.Counter()
        for u in U:
            if not (u.get('lemma') and RS.eligible(u['lemma'])): continue
            for r in u['readings']:
                if RS.eligible(r): cm[f"lemma={RS.spelling_vs(r, u['lemma'], broad)} window={L.spelling(u, r, broad)}"] += 1
        ag['broad' if broad else 'strict'] = dict(cm)
    out['agreement'][t] = ag
    sh = [(u, j) for u in U for j, g in enumerate(u['groups'][1:]) if len(set(g) & Ws) >= 2]
    out['lemma_coverage'][t] = dict(shared_departures=len(sh), with_printed_lemma=sum(1 for u, _ in sh if u.get('lemma')))
    for broad in (False, True):
        for u, j in sh:
            r = u['readings'][j]
            if RS.classify_unit(u, L, broad)[j]:
                if u.get('lemma'): m = u['lemma']
                else:
                    W = [x.strip(',.;:!?|()[]"\'') for x in RS.words(re.sub(r'(\w)-\s+(\w)', r'\1\2', ' '.join(L.window(u, r))))]
                    na = [RS.norm(x, broad) for x in RS.clean_reading(r)]; k = len(na)
                    m = next((' '.join(W[i:i + k]) for i in range(len(W) - k + 1) if [RS.norm(x, broad) for x in W[i:i + k]] == na), '?')
                samp['broad' if broad else 'strict'].append(dict(text=t, leaf=u['leaf'], line=u['entry_line'], group=''.join(sorted(set(u['groups'][j + 1]) & Ws)), reading=r, matched=m, lemma_printed=bool(u.get('lemma'))))
rng = random.Random(20261123)
for k, v in samp.items(): out['sample'][k] = dict(n=len(v), drawn=rng.sample(v, min(40, len(v))))
json.dump(out, open(os.path.join(ROOT, 'checks_20261123', 'spelling_check.json'), 'w'), ensure_ascii=False, indent=1)
tot = collections.defaultdict(collections.Counter)
for t, a in out['agreement'].items():
    for k, c in a.items(): tot[k].update(c)
print({k: dict(v) for k, v in tot.items()})
for k in out['sample']:
    print('==', k, out['sample'][k]['n'])
    for s in out['sample'][k]['drawn']: print(f"   {s['text'][:14]:14s} {s['group']:6s} {s['reading']!r:28s} ~ {s['matched']!r}")
