"""OCR dMGH pages with Transkribus (Metagrapho v1, Text Titan htrId 356425, ~0.5 credit/page).
Usage: tk_ocr.py OUTDIR VOL:LEAF [VOL:LEAF ...]   (or --list file with one VOL:LEAF per line)
Writes OUTDIR/VOL/LEAF.txt (text lines in reading order). Skips pages already done."""
import sys, os, time, json, re, urllib.request, concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import transkribus_api as tk
from xml.etree import ElementTree as ET
out = sys.argv[1]; items = sys.argv[2:]
if items and items[0] == '--list': items = [l.strip() for l in open(items[1]) if l.strip()]
TOK = {'t': tk.token(), 'at': time.time()}
def tok():
    if time.time() - TOK['at'] > 240: TOK['t'] = tk.token(); TOK['at'] = time.time()
    return TOK['t']
def one(it):
    vol, leaf = it.split(':'); leaf = int(leaf)
    fn = os.path.join(out, vol, '%05d.txt' % leaf)
    if os.path.exists(fn) and os.path.getsize(fn) > 0: return it + ' skip'
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    url = 'https://www.dmgh.de/%s/img/300/%s_%05d.jpg' % (vol, vol, leaf)
    r = tk.api(tok(), 'POST', 'processes', {'config': {'textRecognition': {'htrId': 356425}}, 'image': {'imageUrl': url}})
    pid = r['processId']
    for _ in range(120):
        time.sleep(5)
        s = tk.api(tok(), 'GET', 'processes/%s' % pid)
        if s.get('status') in ('FINISHED', 'FAILED', 'CANCELED'): break
    if s.get('status') != 'FINISHED': return '%s %s %s' % (it, pid, s.get('status'))
    req = urllib.request.Request('https://transkribus.eu/processing/v1/processes/%s/page' % pid,
                                 headers={'Authorization': 'Bearer ' + tok()})
    xml = urllib.request.urlopen(req, timeout=120).read()
    root = ET.fromstring(xml); ns = {'p': root.tag.split('}')[0].strip('{')}
    lines = []
    for tl in root.iter('{%s}TextLine' % ns['p']):
        u = tl.find('.//p:TextEquiv/p:Unicode', ns)
        if u is not None and u.text: lines.append(u.text)
    open(fn, 'w').write('\n'.join(lines) + '\n')
    return '%s %s ok %d lines' % (it, pid, len(lines))
with cf.ThreadPoolExecutor(max_workers=int(os.environ.get('TK_WORKERS', '6'))) as ex:
    for line in ex.map(lambda it: (lambda: one(it))() if True else None, items):
        print(line, flush=True)
