"""Assemble per-leaf Tesseract OCR (dMGH 300 dpi) into one text per edition, keeping only the edition's pages
(selected by running head and page range) and marking each page with '@@LEAF n'.
Usage: reocr_assemble.py OCRDIR VOLUME OFFSET HEAD_REGEX PAGE_FROM PAGE_TO OUT.txt
"""
import sys, os, re
d, vol, off, head, p0, p1, out = sys.argv[1:8]; off, p0, p1 = int(off), int(p0), int(p1)
n = 0
with open(out, 'w') as f:
    for leaf in range(p0 + off, p1 + off + 1):
        fn = os.path.join(d, vol, '%05d.txt' % leaf)
        if not os.path.exists(fn): continue
        L = open(fn, encoding='utf-8', errors='replace').read().split('\n')
        top = ' '.join(x for x in L[:3] if x.strip())
        if not re.search(head, top, re.I): continue
        f.write('@@LEAF %d\n' % leaf); f.write('\n'.join(L) + '\n'); n += 1
print(out, n, 'pages')
