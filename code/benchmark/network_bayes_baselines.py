"""Two further baselines, scored exactly like the others (Table I):

1. Neighbor-Net split network (phangorn::neighborNet on pairwise Hamming distances with pairwise deletion of '?').
   The network is not a tree; to score it we take its splits in decreasing weight and keep each one compatible with
   those already kept (greedy compatible subset), which gives a tree.
2. Bayesian inference (MrBayes 3.2.8), Mk model for unordered multistate characters with gamma rate variation,
   2 runs x 4 chains; the majority-rule consensus (posterior >= 0.5) is the reported tree.
   MrBayes allows at most 10 states: states are renumbered per unit by frequency, and states ranked beyond the tenth
   (Heinrichi 91 of 33,108 cells, Donatus 3 of 62,374) are set to missing.

Scoring: Heinrichi = CASC skeleton score against the true tree; Donatus = editor's 12 families recovered exactly.
Usage: network_bayes_baselines.py prep|nnet|bayes|score --ds heinrichi|donatus --out DIR [--ngen N]
"""
import os, sys, json, argparse, subprocess, collections, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'donatus'))
from stemma_mwis import load_num, tree_from_clusters, write_graph

DATA = {'heinrichi': dict(num='data/heinrichi/heinrichi.num', truth='data/heinrichi/correct.txt'),
        'donatus': dict(num='data/donatus/final.num', truth=None)}
RENV = dict(HOME=os.environ['HOME'], PATH='/opt/homebrew/bin:/usr/bin:/bin', LANG='en_US.UTF-8')


def safe(w): return 't_' + re.sub(r'[^A-Za-z0-9]', '_', w)


def prep(ds, out):
    wit, rows = load_num(os.path.expanduser(DATA[ds]['num']))
    os.makedirs(out, exist_ok=True)
    cols = []
    for r in rows:
        c = collections.Counter(x for x in r if x != '?')
        rank = {s: i for i, (s, _) in enumerate(c.most_common())}
        cols.append([str(rank[x]) if x != '?' and rank[x] < 10 else '?' for x in r])
    names = [safe(w) for w in wit]
    json.dump(dict(zip(names, wit)), open(os.path.join(out, 'names.json'), 'w'))
    with open(os.path.join(out, 'matrix.nex'), 'w') as f:
        f.write('#NEXUS\nbegin data;\n dimensions ntax=%d nchar=%d;\n format datatype=standard symbols="0123456789" missing=? gap=-;\n matrix\n' % (len(wit), len(cols)))
        for i, n in enumerate(names): f.write(' %s %s\n' % (n, ''.join(col[i] for col in cols)))
        f.write(' ;\nend;\n')
    with open(os.path.join(out, 'matrix.csv'), 'w') as f:
        f.write('taxon,' + ','.join('c%d' % j for j in range(len(cols))) + '\n')
        for i, n in enumerate(names): f.write(n + ',' + ','.join(col[i] for col in cols) + '\n')


def nnet(ds, out):
    R = r'''
suppressMessages(library(phangorn))
m <- read.csv("%s/matrix.csv", row.names=1, colClasses="character")
x <- lapply(seq_len(nrow(m)), function(i) as.character(unlist(m[i,])))
names(x) <- rownames(m)
d <- phyDat(x, type="USER", levels=as.character(0:9), ambiguity="?")
dm <- as.matrix(dist.hamming(d, ratio=TRUE, exclude="pairwise"))
nmiss <- sum(is.na(dm[upper.tri(dm)]))
dm[is.na(dm)] <- mean(dm[upper.tri(dm)], na.rm=TRUE)   # pairs never observed together: mean distance
nn <- neighborNet(as.dist(dm))
sp <- as.splits(nn)
w <- attr(sp, "weights"); lab <- attr(sp, "labels")
lines <- sapply(seq_along(sp), function(i) paste(format(w[i], digits=10), paste(lab[sp[[i]]], collapse=" "), sep="\t"))
writeLines(c(paste0("#missing_pairs\t", nmiss), lines), "%s/nnet_splits.tsv")
''' % (out, out)
    subprocess.run(['/opt/homebrew/bin/Rscript', '-e', R], check=True, env=RENV)


def bayes(ds, out, ngen):
    blk = ('begin mrbayes;\n set autoclose=yes nowarn=yes seed=2026 swapseed=2026;\n execute matrix.nex;\n'
           ' lset rates=gamma;\n prset statefreqpr=fixed(equal);\n'
           ' mcmc ngen=%d nruns=2 nchains=4 samplefreq=500 printfreq=10000 diagnfreq=10000;\n'
           ' sump; sumt contype=halfcompat conformat=simple;\nend;\n' % ngen)
    open(os.path.join(out, 'run.nex'), 'w').write(blk)
    subprocess.run(['/opt/homebrew/bin/mb', 'run.nex'], cwd=out, check=True,
                   stdout=open(os.path.join(out, 'mb.log'), 'w'), stderr=subprocess.STDOUT)


def clusters_from_newick(s, names, wit):
    """Return the clusters (as sets of witness indices) of a Newick tree, ignoring branch lengths/supports."""
    s = re.sub(r':[0-9.eE+-]+', '', s.strip().rstrip(';'))
    s = re.sub(r'\)[0-9.]+', ')', s)
    out, stack = [], []
    tok = re.findall(r'\(|\)|,|[^(),]+', s)
    cur = [[]]
    for t in tok:
        if t == '(': cur.append([])
        elif t == ')':
            c = cur.pop(); out.append(frozenset(c)); cur[-1].extend(c)
        elif t == ',': pass
        else:
            t = t.strip()
            if t: cur[-1].append(wit.index(names[t]))
    return out


def score(ds, out):
    import build_matrix as bm
    from score import score as casc
    wit, rows = load_num(os.path.expanduser(DATA[ds]['num'])); n = len(wit); FULL = frozenset(range(n))
    names = json.load(open(os.path.join(out, 'names.json')))
    res = {}
    canon = lambda C: FULL - C if 0 in C else C
    # Neighbor-Net -> greedy compatible tree
    if os.path.exists(os.path.join(out, 'nnet_splits.tsv')):
        sp = []
        for ln in open(os.path.join(out, 'nnet_splits.tsv')):
            if ln.startswith('#') or not ln.strip(): continue
            w, mem = ln.rstrip('\n').split('\t')
            sp.append(dict(weight=float(w), members=mem.split()))
        sp = sorted(sp, key=lambda d: -d['weight'])
        kept = []
        for d in sp:
            C = canon(frozenset(wit.index(names[m]) for m in d['members']))
            if not (2 <= len(C) <= n - 2) or C in kept: continue
            if all(not (C & D) or C <= D or D <= C for D in kept): kept.append(C)
        allsplits = [canon(frozenset(wit.index(names[m]) for m in d['members'])) for d in sp]
        res['neighbornet'] = evaluate(ds, kept, allsplits, wit, n, out, 'nnet', bm, casc)
    # MrBayes consensus
    con = os.path.join(out, 'matrix.nex.con.tre')
    if os.path.exists(con):
        t = open(con).read()
        tr = re.search(r'tree\s+con_50_majrule\s*=\s*(?:\[[^\]]*\]\s*)?(.+?);', t, re.S).group(1)
        trans = dict(re.findall(r'^\s*(\d+)\s+(t_\S+?)[,;]?\s*$', t, re.M))
        if trans: tr = re.sub(r'(?<=[(,])\s*(\d+)\s*(?=[,):])', lambda m: trans[m.group(1)], tr)
        cl = [canon(C) for C in clusters_from_newick(tr, names, wit)]
        kept = [C for C in dict.fromkeys(cl) if 2 <= len(C) <= n - 2]
        res['mrbayes'] = evaluate(ds, kept, kept, wit, n, out, 'mrbayes', bm, casc)
        lg = open(os.path.join(out, 'mb.log')).read()
        m = re.findall(r'Average standard deviation of split frequencies:\s*([0-9.]+)', lg)
        res['mrbayes']['asdsf_final'] = float(m[-1]) if m else None
    json.dump(res, open(os.path.join(out, 'scores.json'), 'w'), indent=1)
    print(json.dumps(res, indent=1))


def evaluate(ds, kept, allsplits, wit, n, out, tag, bm, casc):
    r = dict(families_in_tree=len(kept))
    if ds == 'heinrichi':
        g = os.path.join(out, 'tree_%s.txt' % tag)
        write_graph(tree_from_clusters(kept, n), wit, g, directed=False)
        r['skeleton'] = round(casc(os.path.expanduser(DATA[ds]['truth']), g)[0], 2)
    else:
        FULL = frozenset(range(n))
        L = {gn: frozenset(wit.index(w) for w in ms) for gn, ms in bm.GROUPS.items()}
        canonL = {g: (FULL - s if 0 in s else s) for g, s in L.items()}
        r['editor_families_exact'] = sum(1 for s in canonL.values() if s in kept)
        r['which'] = sorted(g for g, s in canonL.items() if s in kept)
        r['editor_families_among_all_splits'] = sum(1 for s in canonL.values() if s in allsplits)
    return r


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('step'); ap.add_argument('--ds', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--ngen', type=int, default=2000000)
    a = ap.parse_args()
    {'prep': lambda: prep(a.ds, a.out), 'nnet': lambda: nnet(a.ds, a.out),
     'bayes': lambda: bayes(a.ds, a.out, a.ngen), 'score': lambda: score(a.ds, a.out)}[a.step]()
