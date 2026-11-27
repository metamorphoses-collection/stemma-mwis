"""Benecke & Lachmann, Iwein (Berlin 1827): apparatus entries -> verified readings with sigla.
Witness sigla (preface): A B C D (parchment), a b c d (paper; d Vienna 1517), e (supplement to B's gaps).
OCR confusions accepted, and only these: J/Ä for A, h for b (inside a run of sigla letters).
A record is kept only if its reading occurs in the entry (normalised) and each siglum occurs in a sigla run of the entry.
Usage: iwein_extract.py iwein_entries.json out.json [--model anthropic/claude-sonnet-5.5]"""
import sys, os, json, re, argparse, urllib.request, concurrent.futures as cf
SIG = set('ABCDabcde')
PROMPT = """Below are apparatus entries from the 1827 edition of Hartmann von Aue's Iwein by Benecke and Lachmann (OCR, with errors). Each entry starts with a line number, then a lemma (sometimes followed by ']'), then variant readings, each followed by the sigla of the manuscripts that have it, separated by ':' or ';'. Manuscript sigla are single letters: A B C D a b c d e. OCR often prints A as 'J' or 'Ä' and b as 'h': read them as A and b. 'fehlt X' means manuscript(s) X omit the lemma: give the reading as "om.". When the lemma itself is followed by sigla (e.g. 'dan ADb: danne diu Bacd'), the lemma is a reading of those sigla.

For each entry return the readings explicitly attributed to sigla, copying each reading exactly as printed (Middle High German words only). Do not infer readings for manuscripts that are not named.

Entries:
{entries}

Answer with JSON only: {{"entries": [{{"verse": 1550, "readings": [{{"reading": "dan", "sigla": ["A","D","b"]}}, {{"reading": "danne diu", "sigla": ["B","a","c","d"]}}]}}]}}"""
norm = lambda t: re.sub(r'[^a-zäöüæœ]', '', t.lower())
def sigla_runs(txt):
    runs = re.findall(r'(?<![A-Za-zäöüß])([ABCDabcdeJÄh]{1,7})(?![A-Za-zäöüß])', txt)
    out = set()
    for r in runs:
        for ch in r: out.add({'J': 'A', 'Ä': 'A', 'h': 'b'}.get(ch, ch))
    return out
def call(model, prompt, key):
    body = dict(model=model, messages=[dict(role='user', content=prompt)], max_tokens=6000, temperature=0, usage=dict(include=True))
    req = urllib.request.Request('https://openrouter.ai/api/v1/chat/completions', data=json.dumps(body).encode(),
                                 headers={'Authorization': f'Bearer {key}', 'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=300) as r: return json.load(r)
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('entries'); ap.add_argument('out'); ap.add_argument('--model', default='anthropic/claude-sonnet-5.5')
    a = ap.parse_args(); key = os.environ['OPENROUTER_API_KEY']
    E = json.load(open(a.entries)); batches = [E[i:i + 25] for i in range(0, len(E), 25)]
    def go(batch):
        txt = '\n'.join(f"{e['verse']}. {e['text']}" for e in batch)
        try:
            o = call(a.model, PROMPT.format(entries=txt), key); t = o['choices'][0]['message'].get('content') or ''
            m = re.search(r'\{.*\}', t, re.S); return batch, json.loads(m.group(0))['entries'], o.get('usage', {}).get('cost', 0)
        except Exception as ex: return batch, [], 0
    out, kept, rej, cost = [], 0, 0, 0
    with cf.ThreadPoolExecutor(8) as ex:
        for batch, res, c in ex.map(go, batches):
            cost += c or 0; byv = {}
            for e in batch: byv.setdefault(e['verse'], []).append(e['text'])
            for r in res:
                texts = byv.get(r.get('verse'))
                if not texts: continue
                full = ' '.join(texts); sr = sigla_runs(full); good = []
                for rd in r.get('readings', []):
                    sg = [s for s in rd.get('sigla', []) if s in SIG]
                    core = rd.get('reading', ''); core = '' if core.startswith('om') else core
                    ok_txt = (core == '' and 'fehlt' in full) or (norm(core) and norm(core) in norm(full))
                    if sg and all(s in sr for s in sg) and ok_txt: good.append(dict(reading=rd['reading'], sigla=sg)); kept += 1
                    else: rej += 1
                if good: out.append(dict(verse=r['verse'], readings=good))
    json.dump(out, open(a.out, 'w'), indent=1, ensure_ascii=False)
    print(f'entries {len(E)}  units with verified readings {len(out)}  readings kept {kept}  rejected {rej}  cost ${cost:.2f}')
