"""Benchmark generative language models on stemma reconstruction (OpenRouter).

Protocol (paper Sec. "Benchmarking against generative models"):
  * identical input for every model: the collation matrix, witnesses anonymised (W01.. in shuffled order,
    seeded) so that no published stemma can be recalled; readings coded per unit;
  * identical task: return nested families (groups descending from one lost exemplar) as JSON;
  * identical scoring: families -> unrooted splits -> tree (Buneman) -> CASC skeleton score (Heinrichi)
    or exact recovery of the editor's 12 families (Donatus); validity = families pairwise compatible
    (crossing families are reported, and dropped greedily in the order given before scoring).
Key: ~/.config/openrouter/env (OPENROUTER_API_KEY), read from the environment, never printed.
Usage: llm_bench.py run --models a,b --datasets heinrichi,donatus --out DIR
       llm_bench.py score --out DIR
"""
import sys, os, json, random, re, time, argparse, urllib.request, concurrent.futures as cf
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'donatus'))
from stemma_mwis import load_num, tree_from_clusters, write_graph

DATA = {
    'heinrichi': dict(num='data/heinrichi/heinrichi.num', truth='data/heinrichi/correct.txt'),
    'donatus': dict(num='data/donatus/final.num', truth=None),
}
CODES = '0123456789abcdefghijklmnopqrstuvwxyz'

PROMPT = """You are an expert in textual criticism (stemmatics).

Below is the full collation of {n} surviving manuscripts of one text. The manuscripts are anonymised as {first}..{last}. Each line gives one manuscript; each character is one variation unit (a place where the manuscripts disagree), in the same order for every line. Within a unit, manuscripts with the same symbol have the same reading; symbols are arbitrary labels and '0' is NOT necessarily the original reading. '?' means the manuscript is not legible or lacks that passage. Units are grouped in blocks of 10 for readability.

{matrix}

Task: reconstruct the stemma codicum, i.e. the genealogy of these manuscripts. Report it as the list of FAMILIES: each family is a set of at least two manuscripts that descend from a common lost exemplar that no manuscript outside the family descends from. Families must be properly nested or disjoint (no two families may cross), as in a tree. Include every family you are reasonably confident of, from large branches down to pairs; do not list the set of all manuscripts.

Answer with your reasoning first if you wish, then end with a single JSON object on its own, exactly in this form:
{{"families": [["W01", "W07"], ["W01", "W07", "W12"], ...]}}"""


def anonymise(wit, seed):
    idx = list(range(len(wit))); random.Random(seed).shuffle(idx)
    width = 2
    alias = {wit[i]: f'W{k + 1:0{width}d}' for k, i in enumerate(idx)}
    return alias


def render(wit, rows, alias):
    inv = sorted(wit, key=lambda w: alias[w])
    lines = []
    for w in inv:
        i = wit.index(w)
        s = ''.join('?' if r[i] == '?' else CODES[int(r[i])] for r in rows)
        lines.append(f'{alias[w]}: ' + ' '.join(s[k:k + 10] for k in range(0, len(s), 10)))
    return '\n'.join(lines)


def call(model, prompt, key, max_tokens=128000):
    # every model at its own 'high' reasoning setting, generous output ceiling; tokens and cost are reported
    body = dict(model=model, messages=[dict(role='user', content=prompt)], max_tokens=max_tokens,
                reasoning=dict(effort='high'), usage=dict(include=True))
    req = urllib.request.Request('https://openrouter.ai/api/v1/chat/completions', data=json.dumps(body).encode(),
                                 headers={'Authorization': f'Bearer {key}', 'Content-Type': 'application/json',
                                          'X-Title': 'stemma-benchmark'})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=3600) as r:
        out = json.load(r)
    out['_seconds'] = round(time.time() - t0, 1)
    return out


def parse_families(text):
    blocks = re.findall(r'\{[^{}]*"families"\s*:\s*\[.*?\]\s*\}', text, flags=re.S)
    for b in reversed(blocks):
        try:
            return json.loads(b)['families']
        except Exception:
            continue
    return None


def run(a):
    key = os.environ['OPENROUTER_API_KEY']
    os.makedirs(a.out, exist_ok=True)
    jobs = []
    for ds in a.datasets.split(','):
        wit, rows = load_num(os.path.expanduser(DATA[ds]['num']))
        alias = anonymise(wit, seed=2026)
        json.dump(alias, open(os.path.join(a.out, f'{ds}_alias.json'), 'w'), indent=0)
        al = sorted(alias.values())
        prompt = PROMPT.format(n=len(wit), first=al[0], last=al[-1], matrix=render(wit, rows, alias))
        open(os.path.join(a.out, f'{ds}_prompt.txt'), 'w').write(prompt)
        for m in a.models.split(','):
            f = os.path.join(a.out, f'{ds}__{m.replace("/", "_")}.json')
            if not os.path.exists(f): jobs.append((ds, m, prompt, f))

    def go(job):
        ds, m, prompt, f = job
        try:
            try:
                out = call(m, prompt, key)
            except Exception as e:                       # model with a smaller output limit: retry once
                if getattr(e, 'code', None) != 400: raise
                out = call(m, prompt, key, max_tokens=32000)
            json.dump(out, open(f, 'w'))
            u = out.get('usage', {})
            return f'{ds:9s} {m:40s} ok  {out["_seconds"]:6.0f}s  cost ${u.get("cost", 0):.3f}'
        except Exception as e:
            msg = str(e)
            if hasattr(e, 'read'):
                try: msg += ' ' + e.read().decode()[:300]
                except Exception: pass
            return f'{ds:9s} {m:40s} FAIL {msg[:300]}'
    with cf.ThreadPoolExecutor(max_workers=a.workers) as ex:
        for line in ex.map(go, jobs): print(line, flush=True)


def score(a):
    from score import score as casc
    from baselines import clusters as _c
    import build_matrix as bm
    res = []
    for f in sorted(os.listdir(a.out)):
        if '__' not in f or not f.endswith('.json'): continue
        ds, m = f[:-5].split('__')
        out = json.load(open(os.path.join(a.out, f)))
        wit, rows = load_num(os.path.expanduser(DATA[ds]['num'])); n = len(wit)
        alias = json.load(open(os.path.join(a.out, f'{ds}_alias.json'))); back = {v: k for k, v in alias.items()}
        msg = out['choices'][0]['message']
        text = msg.get('content') or ''
        fams = parse_families(text)
        u = out.get('usage', {})
        row = dict(dataset=ds, model=m, seconds=out.get('_seconds'), cost=u.get('cost'),
                   tokens_out=u.get('completion_tokens'))
        if fams is None:
            row.update(valid=False, note='no parsable JSON'); res.append(row); continue
        splits, crossing, unknown = [], 0, 0
        FULL = frozenset(range(n))
        for fam in fams:
            idx = frozenset(wit.index(back[x]) for x in fam if x in back)
            unknown += sum(1 for x in fam if x not in back)
            C = FULL - idx if 0 in idx else idx
            if not (2 <= len(C) <= n - 2) or C in splits: continue
            if all(not (C & D) or C <= D or D <= C for D in splits): splits.append(C)
            else: crossing += 1
        row.update(families_given=len(fams), families_kept=len(splits), crossing_dropped=crossing,
                   unknown_ids=unknown, valid=crossing == 0)
        if ds == 'heinrichi':
            g = os.path.join(a.out, f'tree_{ds}__{m}.txt')
            write_graph(tree_from_clusters(splits, n), wit, g, directed=False)
            row['skeleton'] = round(casc(os.path.expanduser(DATA[ds]['truth']), g)[0], 2)
        else:
            L = {gname: frozenset(wit.index(w) for w in ms) for gname, ms in bm.GROUPS.items()}
            canonL = {g: (FULL - s if 0 in s else s) for g, s in L.items()}
            row['editor_families_exact'] = sum(1 for s in canonL.values() if s in splits)
            row['which'] = sorted(g for g, s in canonL.items() if s in splits)
        res.append(row)
    json.dump(res, open(os.path.join(a.out, 'scores.json'), 'w'), indent=1)
    for r in res: print(json.dumps(r))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('--models', default='')
    ap.add_argument('--datasets', default='heinrichi,donatus'); ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    run(a) if a.cmd == 'run' else score(a)
