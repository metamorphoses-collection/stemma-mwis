"""CHECK (2026-11-23): the sample flow 71 -> 46 -> 42 -> 28 -> 21 -> 14, with the count for each
exclusion reason, recomputed from the frozen records (tops/*.json, settings.json, runs/, screen_a10/, proof scores).
Read-only. Output: data/bedier/checks_20261123/flow.json"""
import json, glob, os, collections
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'bedier')
OUT = os.path.join(R, 'checks_20261123'); os.makedirs(OUT, exist_ok=True)
recs = []
for f in sorted(glob.glob(os.path.join(R, 'tops', '*.json'))):
    for t in json.load(open(f))['texts']: recs.append(t)
S = {s['text_id']: s for s in json.load(open(os.path.join(R, 'settings.json')))}
runs = set(os.listdir(os.path.join(R, 'runs')))
fresh = json.load(open(os.path.join(R, 'proof_scores_a10_fresh.json')))
new = json.load(open(os.path.join(R, 'proof_scores_a10_new.json')))
resc = json.load(open(os.path.join(R, 'proof_scores_a10_rescore.json')))
flow, log = collections.OrderedDict(), []
def step(name, keep, reasons):
    flow[name] = dict(n=len(keep), excluded=dict(collections.Counter(r for _, r in reasons)))
    log.extend(dict(stage=name, text=t, reason=r) for t, r in reasons)
flow['recorded'] = dict(n=len(recs)); A = [t['text_id'] for t in recs]
k1 = [t for t in A if t in S]; step('stemma_stated', k1, [(t, 'no stemma (A3 none)') for t in A if t not in S])
k2 = [t for t in k1 if len(S[t]['sigla']) >= 2]
step('two_sigla_and_top', k2, [(t, 'fewer than 2 single-capital sigla (A4 siglum form)') for t in k1 if t not in k2])
k3 = [t for t in k2 if t in runs]; step('located', k3, [(t, 'not located (A7/A8, start or end match < 0.75)') for t in k2 if t not in k3])
k4, r4 = [], []
for t in k3:
    sc = json.load(open(os.path.join(R, 'screen_a10', t + '.json')))
    if sc['included']: k4.append(t)
    else: r4.append((t, ' + '.join(sc['exclusion'])))
step('screen', k4, r4)
k5, r5 = [], []
for t in k4:
    p = fresh.get(t) or new.get(t) or resc.get(t)
    src = 'fresh (A10 repair)' if t in fresh else 'new (A10)' if t in new else 'rescore (A10, parser v4)'
    if p and p['passes']: k5.append(t)
    else: r5.append((t, f'proof page fails [{src}]' + (f" P={p['precision']:.2f} R={p['recall']:.2f}" if p else ' (no score)')))
step('proof_page', k5, r5)
k6, r6 = [], []
for t in k5:
    s = S[t]
    if s['kind'] != 'explicit': r6.append((t, f"top {s['kind']} (reported separately)"))
    elif s['pilot']: r6.append((t, 'pilot / case-study edition'))
    else: k6.append(t)
step('primary', k6, r6)
kinds = collections.Counter(S[t]['kind'] for t in k2)
json.dump(dict(flow=flow, kinds_at_46=kinds, primary=k6, log=log), open(os.path.join(OUT, 'flow.json'), 'w'), indent=1)
for k, v in flow.items(): print(k, v)
print('kinds among the 46:', dict(kinds)); print('primary:', k6)
for l in log: print(f"  {l['stage']:18s} {l['text']:36s} {l['reason']}")
