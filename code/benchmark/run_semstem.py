"""Run Semstem (Roos and Zou 2011) from an unmodified checkout of the Stemweb repository (github.com/Stemweb/Stemweb,
Stemweb/algorithms/semstem/semstem.py), outside the Stemweb web framework. The Python 2 source is loaded and patched
mechanically at run time (print statements, dict views, framework stubs); the algorithm is not changed. Stemweb's code
carries no licence, so it is not copied here.

Input: the NEXUS matrix written by network_bayes_baselines.py prep (states 0-9, '?' missing).
Output: OUT/semstem_seed<k>/tree/treebest.tre, then scored by network_bayes_baselines.py score.
usage: run_semstem.py STEMWEB_CHECKOUT --nex DIR/matrix.nex --out DIR --seeds 1,2,3 --iter 200
"""
import os, re, sys, types, random, argparse
import numpy as np


def load_semstem(checkout):
    src = open(os.path.join(checkout, 'Stemweb', 'algorithms', 'semstem', 'semstem.py')).read()
    src = src.replace('from Stemweb.algorithms.tasks import AlgorithmTask', 'class AlgorithmTask(object): pass')
    src = src.replace('from Stemweb.algorithms.utils import newick2img', '')
    src = src.replace('from numpy import *', 'from numpy import *\nInf = inf; NaN = nan\nimport builtins as _b; max = _b.max; min = _b.min', 1)   # NumPy 2 removed these aliases
    src = re.sub(r'^(\s*)utils\.tree2img\(.*$', r'\1pass', src, flags=re.M)
    src = re.sub(r"^(\s*)os\.system\('neato.*$", r'\1pass', src, flags=re.M)              # drawing only
    src = src.replace('nodelistnotextendedtemp.remove(node)',
                      '(node in nodelistnotextendedtemp and nodelistnotextendedtemp.remove(node))')
    src = src.replace('nodelistnotextendedtemp.remove(neighbori)',                       # export guard (dict order)
                      '(neighbori in nodelistnotextendedtemp and nodelistnotextendedtemp.remove(neighbori))')
    src = re.sub(r'^(\s*)print (.+)$', r'\1print(\2)', src, flags=re.M)
    src = re.sub(r'(\w+)\.has_key\(([^)]*)\)', r'(\2 in \1)', src)
    src = re.sub(r'(\w+)\.keys\(\)', r'list(\1.keys())', src)
    src = re.sub(r'(\w+)\.values\(\)', r'list(\1.values())', src)
    src = re.sub(r'(\w+)\.items\(\)', r'list(\1.items())', src)
    mod = types.ModuleType('semstem_port'); exec(compile(src, 'semstem.py', 'exec'), mod.__dict__)
    return mod.Semstem


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('checkout'); ap.add_argument('--nex', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--seeds', default='1,2,3'); ap.add_argument('--iter', type=int, default=200)
    a = ap.parse_args(); Semstem = load_semstem(a.checkout)
    for sd in [int(x) for x in a.seeds.split(',')]:
        random.seed(sd); np.random.seed(sd)
        out = os.path.abspath(os.path.join(a.out, 'semstem_seed%d' % sd)); os.makedirs(out, exist_ok=True)
        t = Semstem(); t.run_args = dict(infile=os.path.abspath(a.nex), itermaxin=a.iter, outfolder=out, learnlength=False)
        t.newick_path = os.path.join(out, 'tree', 'treebest.tre'); t.image_path = os.path.join(out, 'tree', 'treebest.png')
        t._put_in_results_ = lambda d: open(os.path.join(out, 'qscore.txt'), 'w').write(str(d))
        t._stop = types.SimpleNamespace(value=0)
        sys.stdout = open(os.path.join(out, 'stdout.log'), 'w')
        try: t.__algorithm__(t.run_args)
        finally: sys.stdout.close(); sys.stdout = sys.__stdout__
        std = standardise(open(t.newick_path).read())
        open(os.path.join(out, 'tree', 'treebest_std.tre'), 'w').write(std)
        print('seed', sd, open(os.path.join(out, 'qscore.txt')).read())


def standardise(newick):
    """Semstem writes surviving witnesses as internal nodes, e.g. ((()t_L)t_K)t_B: rewrite so that a labelled internal node
    becomes a leaf beside its children, '(children,LABEL)', and drop the empty '()' leaves."""
    s = newick.strip().rstrip(';').replace('()', '')
    out, i = [], 0
    def node(i):
        if s[i] == '(':
            kids = []; i += 1
            while True:
                k, i = node(i); kids.append(k)
                if s[i] == ',': i += 1; continue
                if s[i] == ')': i += 1; break
            j = i
            while j < len(s) and s[j] not in ',();': j += 1
            lab = s[i:j]; kids = [k for k in kids if k]
            if lab: kids.append(lab)
            return ('(' + ','.join(kids) + ')') if len(kids) > 1 else (kids[0] if kids else ''), j
        j = i
        while j < len(s) and s[j] not in ',();': j += 1
        return s[i:j], j
    t, _ = node(0)
    return t + ';'


if __name__ == '__main__':
    main()
