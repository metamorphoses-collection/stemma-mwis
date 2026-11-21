"""Export the machine-read critical apparatus of the MGH Auctores antiquissimi texts as plain tables.
usage: export_apparatus_tables.py BEDIER_DATA_DIR OUTDIR
BEDIER_DATA_DIR holds runs_a10/<text>/{units.json,leaves.json,llm.txt}, tops/*.json, classification_a10/,
proof_scores*.json. Writes OUTDIR/apparatus/<text>.csv, texts.csv, sigla.csv, proof_pages.csv."""
import json, csv, glob, os, sys
B, O = sys.argv[1], sys.argv[2]
os.makedirs(f'{O}/apparatus', exist_ok=True)
meta = {}
for f in sorted(glob.glob(f'{B}/tops/*.json')):
    for t in json.load(open(f)).get('texts', []): meta[t['text_id']] = t
texts, sig_rows, n_all = [], [], 0
for d in sorted(glob.glob(f'{B}/runs_a10/*/units.json')):
    tid = d.split('/')[-2]; R = os.path.dirname(d)
    U = json.load(open(d)); L = open(f'{R}/llm.txt').read().split('\n'); lv = json.load(open(f'{R}/leaves.json'))
    m = meta.get(tid, {}); c = {}
    cf = f'{B}/classification_a10/{tid}.json'
    if os.path.exists(cf): c = json.load(open(cf))
    n = 0
    with open(f'{O}/apparatus/{tid}.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['entry', 'volume', 'image_leaf', 'text_line', 'lemma', 'reading', 'witnesses', 'apparatus_line'])
        for k, u in enumerate(U, 1):
            line = L[u['entry_line'] - 1].strip() if 0 < u['entry_line'] <= len(L) else ''
            for g, r in zip(u['groups'][1:], u['readings']):
                w.writerow([k, lv['vol'], u['leaf'], u['line'], u.get('lemma') or '', r, g, line]); n += 1
    n_all += n
    texts.append(dict(text_id=tid, work=m.get('work', ''), editor=m.get('editor', ''), year=m.get('year', ''),
        volume=m.get('volume', lv['vol']), dmgh_volume=lv['vol'], first_leaf=lv['leaves'][0], last_leaf=lv['leaves'][1],
        apparatus_entries=len(U), readings=n, witnesses_cited=''.join(sorted(c.get('cited', {}), key=lambda s: -c['cited'][s])) if c else '',
        editor_top_division=' | '.join(c.get('editor_top', [])) if c else ''))
    for s in m.get('sigla', []):
        sig_rows.append(dict(text_id=tid, siglum=s.get('siglum', ''), manuscript=s.get('manuscript', ''),
                             date=s.get('date_if_given', '') or '', role=s.get('role', '') or ''))
def dump(name, rows):
    with open(f'{O}/{name}', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
dump('texts.csv', texts); dump('sigla.csv', sig_rows)
pr = []
for f in sorted(glob.glob(f'{B}/proof_scores*.json')):
    rnd = os.path.basename(f)[len('proof_scores'):-5].strip('_') or 'a1'
    for tid, s in json.load(open(f)).items():
        pr.append(dict(text_id=tid, round=rnd, printed_page=s.get('page', ''), image_leaf=s.get('leaf', ''),
                       groups_on_page=s.get('printed', ''), recovered=s.get('right', ''), extra=s.get('extra', ''),
                       precision=round(s['precision'], 3) if s.get('precision') is not None else '',
                       recall=round(s['recall'], 3) if s.get('recall') is not None else '', passes=s.get('passes', '')))
dump('proof_pages.csv', pr)
print(len(texts), 'texts', n_all, 'readings', len(sig_rows), 'sigla rows', len(pr), 'proof rows')
