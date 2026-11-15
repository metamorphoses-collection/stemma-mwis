"""Roelli and Bachmann's Leitfehler method (Rev. hist. textes 5, 2010, 307-331), as described in the Handbook of
Stemmatology (2020), sect. 5.3.7.2 and 5.5.4, rebuilt here as a baseline; scored like the others (Table I).

1. Every reading of a unit becomes a binary character (has it / has another reading / missing): absence-presence form.
2. Candidates are the informative characters (at least two witnesses with the reading and two without).
3. Two candidates are compatible unless all four combinations (1,1), (1,0), (0,1), (0,0) occur among the witnesses
   observed in both (the pairwise test). Each candidate's count = number of other candidates compatible with it.
4. Score = SMAX * count / max count (SMAX = 20; the tree is invariant to this scale).
5. Distance between two witnesses = score-weighted share of differing candidates over the candidates observed in both
   (the normalisation over co-observed characters is ours: the method was designed for complete witnesses).
6. Tree: Fitch-Margoliash weighted least squares (weights 1/d^2), searched by nearest-neighbour interchanges from the
   neighbour-joining tree (R phangorn nnls.tree; PHYLIP is not used). Not reproduced: the down-weighting of runs of k
   consecutive syntactic edits (parameter p), which needs word-level alignment, and Roelli's later hand-picked stages.

Usage: leitfehler_rb.py --ds heinrichi|donatus --out DIR      (then network_bayes_baselines.py score --ds .. --out DIR)
"""
import os, sys, json, argparse, subprocess, itertools, collections, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from stemma_mwis import load_num
from network_bayes_baselines import DATA, RENV, safe

SMAX = 20.0


def characters(rows, n):
    C = []
    for r in rows:
        states = collections.Counter(x for x in r if x != '?')
        if len(states) < 2: continue
        for s in states:
            v = np.array([-1 if x == '?' else (1 if x == s else 0) for x in r], dtype=np.int8)
            if (v == 1).sum() >= 2 and (v == 0).sum() >= 2: C.append(v)
    # one character per bipartition and missing pattern: a binary unit gives the same split twice
    seen, out = set(), []
    for v in C:
        flip = np.where(v == -1, -1, 1 - v)
        key = min(tuple(v.tolist()), tuple(flip.tolist()))
        if key not in seen: seen.add(key); out.append(v)
    return np.array(out)


def compat_counts(X):
    m = len(X); obs = (X >= 0); one = (X == 1); zero = (X == 0)
    cnt = np.zeros(m, dtype=np.int64)
    O = obs.astype(np.int32); A1 = one.astype(np.int32); A0 = zero.astype(np.int32)
    n11 = A1 @ A1.T; n10 = A1 @ A0.T; n01 = A0 @ A1.T; n00 = A0 @ A0.T
    incompatible = (n11 > 0) & (n10 > 0) & (n01 > 0) & (n00 > 0)
    np.fill_diagonal(incompatible, False)
    return (m - 1) - incompatible.sum(1)


def distances(X, w, n):
    D = np.zeros((n, n)); undefined = 0
    for i, j in itertools.combinations(range(n), 2):
        both = (X[:, i] >= 0) & (X[:, j] >= 0)
        tot = w[both].sum()
        if tot == 0: D[i, j] = D[j, i] = np.nan; undefined += 1; continue
        D[i, j] = D[j, i] = (w[both] * (X[both, i] != X[both, j])).sum() / tot
    if undefined: D[np.isnan(D)] = np.nanmean(D[np.triu_indices(n, 1)])
    return D, undefined


R_FM = r'''
suppressMessages(library(phangorn))
d <- as.matrix(read.csv("%(out)s/rb_dist.csv", row.names=1, check.names=FALSE))
dm <- as.dist(d)
W <- 1 / pmax(d, 1e-3)^2; diag(W) <- 0; wd <- as.vector(as.dist(W))
fit <- function(t) { f <- nnls.tree(dm, t, method="unrooted", weight=wd, trace=0)
  o <- cophenetic(f)[rownames(d), rownames(d)]; list(t=f, ss=sum(wd * (as.vector(as.dist(o)) - as.vector(dm))^2)) }
cur <- fit(nj(dm)); best <- cur$ss; it <- 0
repeat { it <- it + 1; nb <- nni(cur$t); improved <- FALSE
  for (k in seq_along(nb)) { cand <- fit(nb[[k]]); if (cand$ss < best - 1e-12) { best <- cand$ss; cur <- cand; improved <- TRUE } }
  if (!improved || it > 200) break }
write.tree(cur$t, "%(out)s/rb_fm.tre")
writeLines(c(paste0("weighted_ss\t", best), paste0("nni_rounds\t", it)), "%(out)s/rb_fm.log")
'''


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--ds', required=True); ap.add_argument('--out', required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    wit, rows = load_num(os.path.expanduser(DATA[a.ds]['num'])); n = len(wit)
    X = characters(rows, n); cnt = compat_counts(X); w = SMAX * cnt / cnt.max()
    D, und = distances(X, w, n)
    names = [safe(x) for x in wit]
    json.dump(dict(zip(names, wit)), open(os.path.join(a.out, 'names.json'), 'w'))
    with open(os.path.join(a.out, 'rb_dist.csv'), 'w') as f:
        f.write(',' + ','.join(names) + '\n')
        for i in range(n): f.write(names[i] + ',' + ','.join('%.10f' % D[i, j] for j in range(n)) + '\n')
    subprocess.run(['/opt/homebrew/bin/Rscript', '-e', R_FM % dict(out=a.out)], check=True, env=RENV)
    info = dict(witnesses=n, candidates=int(len(X)), max_count=int(cnt.max()), median_count=float(np.median(cnt)),
                share_fully_compatible=float((cnt == cnt.max()).mean()), pairs_never_coobserved=und,
                fm=dict(l.strip().split('\t') for l in open(os.path.join(a.out, 'rb_fm.log'))))
    json.dump(info, open(os.path.join(a.out, 'rb_info.json'), 'w'), indent=1); print(json.dumps(info))


if __name__ == '__main__':
    main()
