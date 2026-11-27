"""Run RHM (Roos, Heikkilä and Myllymäki 2006) from the Stemweb repository (Stemweb/algorithms/rhm/binarysankoff_python.c),
compiled standalone (the Python wrapper is cut, the program's own main() is used) with one patch: the original measured
compressed length through a cast of zlib's gzFile to z_stream, which reads zero with current zlib; a helper computes the
same quantity (raw-deflate output after a sync flush). Stemweb's code carries no licence, so it is not copied here.

Input: the aligned CASC table (heinrichi_aligned.csv, tab-separated, one witness per column, empty cell = no word), written
as one file per witness with one aligned slot per line, which is RHM's input format.
RHM seeds itself from the clock, so runs are repeated; each writes sankoff-tree_0.dot, converted here to Newick.
usage: run_rhm.py build STEMWEB_CHECKOUT BINARY_OUT          (patch and compile)
       run_rhm.py RHM_BINARY ALIGNED_CSV --out DIR --runs 3 --iter 100000
"""
import os, re, sys, json, argparse, subprocess, collections


def write_inputs(csv, d):
    rows = [l.rstrip('\n').split('\t') for l in open(csv, encoding='latin-1')]
    hdr = [h.replace('.txt', '') for h in rows[0]]; os.makedirs(d, exist_ok=True)
    for k, h in enumerate(hdr):
        with open(os.path.join(d, h), 'w', encoding='latin-1') as f:
            for r in rows[1:]: f.write((r[k] if k < len(r) else '') + '\n')
    return hdr


def dot_to_newick(dot):
    lab = {int(a): b for a, b in re.findall(r'^(\d+) \[label="([^":]+)"', dot, re.M)}
    adj = collections.defaultdict(list)
    for a, b in re.findall(r'^(\d+) -- (\d+)', dot, re.M): adj[int(a)].append(int(b)); adj[int(b)].append(int(a))
    internal = [v for v in adj if v not in lab]
    root = internal[0] if internal else next(iter(adj))
    def rec(v, p):
        kids = [rec(w, v) for w in adj[v] if w != p]
        if v in lab: return lab[v] if not kids else '(' + ','.join(kids + [lab[v]]) + ')'
        return '(' + ','.join(kids) + ')'
    return rec(root, None) + ';'


HELPER = r"""
/* Patch: the original read the compressed length through a cast of gzFile to z_stream, which only worked with old zlib
   internals. zsize() computes the same quantity, the raw-deflate output of buf after a sync flush. */
static int zsize(const char *s)
{
  z_stream z; unsigned char out[65536]; memset(&z, 0, sizeof z);
  deflateInit2(&z, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY);
  z.next_in = (Bytef *) s; z.avail_in = (uInt) strlen(s);
  do { z.next_out = out; z.avail_out = sizeof out; deflate(&z, Z_SYNC_FLUSH); } while (z.avail_out == 0);
  int t = (int) z.total_out; deflateEnd(&z); return t;
}
"""


def build(checkout, binary):
    src = open(os.path.join(checkout, 'Stemweb', 'algorithms', 'rhm', 'binarysankoff_python.c')).read()
    src = src[:src.index('Simple wrapper for using this module from python')].rsplit('/*', 1)[0]
    src = src.replace('#include <Python.h>', '')
    assert src.count('(int) ((struct z_stream_s *)gfile)->total_out') == 3
    src = src.replace('(int) ((struct z_stream_s *)gfile)->total_out', 'zsize(buf)')
    i = src.index('const int GZIP_HEADER'); j = src.index('\n', i) + 1
    src = src[:j] + HELPER + src[j:]
    c = binary + '.c'; open(c, 'w').write(src)
    subprocess.run(['cc', '-O2', '-w', '-o', binary, c, '-lz'], check=True)


def main():
    if sys.argv[1:2] == ['build']:
        build(sys.argv[2], sys.argv[3]); return
    ap = argparse.ArgumentParser(); ap.add_argument('binary'); ap.add_argument('csv'); ap.add_argument('--out', required=True)
    ap.add_argument('--runs', type=int, default=3); ap.add_argument('--iter', type=int, default=100000); a = ap.parse_args()
    out = os.path.abspath(a.out); inp = os.path.join(out, 'rhm_input'); hdr = write_inputs(a.csv, inp)
    json.dump({h: h for h in hdr}, open(os.path.join(out, 'names.json'), 'w'))
    for r in range(1, a.runs + 1):
        d = os.path.join(out, 'rhm_run%d' % r); os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'stdout.log'), 'w') as so, open(os.path.join(d, 'stderr.log'), 'w') as se:
            subprocess.run([os.path.abspath(a.binary), inp, '.', str(a.iter), '1'], cwd=d, stdout=so, stderr=se, check=True)
        dot = open(os.path.join(d, 'sankoff-tree_0.dot')).read()
        score = re.search(r'sankoff-score (\d+)', dot).group(1)
        open(os.path.join(d, 'rhm.tre'), 'w').write(dot_to_newick(dot))
        open(os.path.join(d, 'score.txt'), 'w').write(score)
        print('run', r, 'sankoff-score', score)


if __name__ == '__main__':
    main()
