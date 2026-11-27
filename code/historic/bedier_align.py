"""A5: locate each proof page in the Internet Archive OCR by its content, not by page numbers.
For the proof leaf and the next leaf, Apple Vision reads the page image; the first lines of edited text (below the page
head, at least 4 words) are fuzzy-matched (difflib ratio on letters only) against the OCR lines of the text's range.
Window = [line matching this page's first text lines, line matching the next page's first text lines), the latter searched
within 30-250 lines after the start; located only if both matches have ratio >= 0.75 (else the proof fails, conservatively).
Independent of the gold transcription and of the parser.
Usage: bedier_align.py LEAVES.json VISION_DIR --out windows.json
"""
import json, os, re, sys, difflib

B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'bedier')
norm = lambda s: re.sub(r'[^a-z]', '', s.lower())


def first_text_lines(js, k=4):
    o = json.load(open(js))['observations']
    H = json.load(open(js)).get('height', 3508)
    body = [x for x in sorted(o, key=lambda x: (round(x['y0'] / 25), x['x0'])) if x['y0'] > 0.12 * H and len(x['text'].split()) >= 4]
    return [x['text'] for x in body[:k]]


def locate(lines, lo, hi, anchors, after=0, until=None):
    best = (0, None)
    N = [norm(l) for l in lines]
    for rank, a in enumerate(anchors):
        na = norm(a)
        if len(na) < 12: continue
        for i in range(max(lo - 1, after), min(hi, len(lines), until or 10**9)):
            if not N[i] or abs(len(N[i]) - len(na)) > 0.6 * len(na): continue
            r = difflib.SequenceMatcher(None, na, N[i]).ratio()
            if r > best[0]: best = (r, i - rank)
    return best


def main():
    leaves, vdir, out = sys.argv[1], sys.argv[2], sys.argv[4]
    st = {s['text_id']: s for s in json.load(open(os.path.join(B, 'settings.json')))}
    R = json.load(open(leaves)); W = {}
    for t, r in R.items():
        s = st[t]; lines = open(os.path.join(B, 'ocr', s['ia_id'] + '.txt'), encoding='utf-8', errors='replace').readlines()
        lo, hi = s['lines']
        a = first_text_lines(os.path.join(vdir, r['vol'], f"{r['leaf']:05d}.v.json"))
        b = first_text_lines(os.path.join(vdir, r['vol'], f"{r['leaf'] + 1:05d}.v.json"))
        ra, ia = locate(lines, lo, hi, a)
        rb, ib = locate(lines, lo, hi + 400, b, after=(ia or 0) + 30, until=(ia or 0) + 250)   # a page spans 30-250 OCR lines
        ok = ra >= 0.75 and rb >= 0.75 and ia is not None and ib is not None
        W[t] = dict(page=r['page'], located=ok, start=None if not ok else ia + 1, end=None if not ok else ib + 1,
                    ratio_start=round(ra, 2), ratio_end=round(rb, 2), anchors_start=a[:2], anchors_end=b[:2])
        print(f"{t[:32]:32s} p.{r['page']:>4} lines {W[t]['start']}-{W[t]['end']}  match {ra:.2f}/{rb:.2f}")
    json.dump(W, open(out, 'w'), ensure_ascii=False, indent=1)


if __name__ == '__main__':
    main()
