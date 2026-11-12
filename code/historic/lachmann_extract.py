"""Lachmann, Lucretius commentary (Berlin 1850): turn verse notes into structured readings, with verification.

Each note (verse number + Latin prose) is sent to a model, which returns, for every reading the note attributes
to a source, {"reading": ..., "sources": [...]}. Sources are normalised to:
  O = oblongus, O2 = corrector oblongi, Q = quadratus, Q2 = corrector quadrati, S = schedae (Gottorpienses/Vindobonenses),
  It = Itali / codices Itali, ED = early printed editions (Brixiensis, Veronensis, Veneta, Aldina, Iuntina ...),
  plus named scholars (Marullus, Lambinus, ...) kept by name but not used as witnesses.
VERIFICATION (no hallucinated evidence): a record is kept only if (i) every source keyword occurs in the note
and (ii) the reading, after normalisation (lowercase, letters only), occurs in the note. Rejected records are logged.
Usage: lachmann_extract.py notes.json out.json --book 1 [--model anthropic/claude-sonnet-5.5]
"""
import sys, os, json, re, time, argparse, urllib.request, concurrent.futures as cf
KEYS = {'O': r'oblong', 'O2': r'corrector\w* oblong|oblong\w* corrector', 'Q': r'quadrat', 'Q2': r'corrector\w* quadrat|quadrat\w* corrector',
        'S': r'sched', 'It': r'ital', 'ED': r'editio|edition|impress|venet|brix|veron|ald|iunt|avant'}
PROMPT = """You extract manuscript evidence from Karl Lachmann's Latin commentary on Lucretius (1850). For the note below, list every reading that the note explicitly attributes to a source. Sources: 'oblongus' -> O; 'corrector oblongi' -> O2; 'quadratus' -> Q; 'corrector quadrati' -> Q2; 'schedae' (Gottorpienses, Vindobonenses) -> S; 'Itali' / Italian manuscripts -> It; early printed editions (Brixiensis, Veronensis, Veneta, Aldina, Iuntina, 'olim impressi', 'editiones') -> ED; named scholars (Marullus, Lambinus, Turnebus, Pius, Avancius, Gifanius, Lachmann's own proposal) -> their name. Copy each reading EXACTLY as printed in the note (a word or phrase, Latin only, no comments). If the note says a source omits a word, the reading is "om. <word>". Do not infer anything that is not stated. If the note attributes no reading to any source, return an empty list.

Note on verse {verse}:
{text}

Answer with JSON only: {{"readings": [{{"reading": "...", "sources": ["O", "Q", ...]}}]}}"""

def call(model, prompt, key):
    body = dict(model=model, messages=[dict(role='user', content=prompt)], max_tokens=1500, temperature=0, usage=dict(include=True))
    req = urllib.request.Request('https://openrouter.ai/api/v1/chat/completions', data=json.dumps(body).encode(),
                                 headers={'Authorization': f'Bearer {key}', 'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=300) as r: return json.load(r)

norm = lambda t: re.sub(r'[^a-z]', '', t.lower().replace('j', 'i').replace('v', 'u'))

def verify(note, rec):
    txt = note['text']; ok_src = []
    for s in rec.get('sources', []):
        k = KEYS.get(s)
        if k is None: continue                      # scholars: kept out of the witness matrix
        if re.search(k, txt, re.I): ok_src.append(s)
    rd = rec.get('reading', '')
    core = rd[4:] if rd.startswith('om. ') else rd
    found = bool(norm(core)) and norm(core) in norm(txt)
    return ok_src, found

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('notes'); ap.add_argument('out'); ap.add_argument('--book', type=int, default=1)
    ap.add_argument('--model', default='anthropic/claude-sonnet-5.5'); ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args(); key = os.environ['OPENROUTER_API_KEY']
    notes = [n for n in json.load(open(a.notes)) if n['book'] == a.book]
    def go(n):
        try:
            o = call(a.model, PROMPT.format(verse=n['verse'], text=n['text'][:6000]), key)
            t = o['choices'][0]['message'].get('content') or ''
            m = re.search(r'\{.*\}', t, re.S); recs = json.loads(m.group(0))['readings'] if m else []
            return n, recs, o.get('usage', {}).get('cost', 0)
        except Exception as e:
            return n, None, 0
    out, cost, kept, rej = [], 0, 0, 0
    with cf.ThreadPoolExecutor(a.workers) as ex:
        for n, recs, c in ex.map(go, notes):
            cost += c or 0
            if recs is None: out.append(dict(verse=n['verse'], error=True)); continue
            good, bad = [], []
            for r in recs:
                src, found = verify(n, r)
                if found and src: good.append(dict(reading=r['reading'], sources=src)); kept += 1
                elif src: bad.append(r); rej += 1
            out.append(dict(verse=n['verse'], line=n['line'], readings=good, rejected=bad))
    json.dump(out, open(a.out, 'w'), indent=1, ensure_ascii=False)
    print(f'notes {len(notes)}  witness readings kept {kept}  rejected {rej}  cost ${cost:.2f}')
