"""MWIS protocol v2: exact MWIS backbone + bootstrap stability + parsimony refinement under the backbone.

1. Exact partial-data MWIS (exact_partial.solve) on the full data.
2. B bootstrap replicates over variation units -> stability of every family (clusters oriented away from witness 0).
3. For each threshold theta: backbone = families with stability >= theta (greedy by stability, kept only if
   compatible with those already kept).
4. Parsimony (phangorn pratchet) on the data plus K heavy constraint characters per backbone family
   (family = 1, rest = 0): parsimony may not break a backbone family, but resolves everything else.
Outputs per theta: backbone size, families of the final tree, score (Heinrichi: CASC skeleton vs truth;
or any scorer supplied by the caller), and a stability table of all families seen (>= 10%).
Usage: mwis_protocol.py <matrix.num> --out DIR [--B 50] [--thetas 0.5,0.7,0.9] [--truth correct.txt]
"""
import sys, os, json, random, subprocess, argparse, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'donatus'))
from stemma_mwis import load_num, tree_from_clusters, write_graph
import exact_partial
from baselines import parse_newick, clusters as tree_clusters

RENV = dict(HOME=os.environ['HOME'], PATH=os.environ.get('PATH', '/usr/bin:/bin'), LANG='en_US.UTF-8')


def canon(C, n):
    C = frozenset(C)
    if 0 in C: C = frozenset(range(n)) - C
    return C if 2 <= len(C) <= n - 2 else None


REPLICATES = []   # per solve: status, objective, bound, time limit used (saved in protocol.json)


def mwis_families(wit, rows, tl, tag='full', tries=3):
    """Families of the certified exact optimum. On a timeout the solve is repeated with four times the limit;
    a replicate that is still uncertified is recorded and contributes no families (never an uncertified tree)."""
    n = len(wit); t = tl
    for k in range(tries):
        try:
            summ, cl, _ = exact_partial.solve(wit, rows, t, log=lambda *a: None)
            REPLICATES.append(dict(tag=tag, status=summ['master_status'], obj=summ['exact_obj'], bound=summ['bound'],
                                   tl=t, certified=True))
            print(f"{tag}: {summ['master_status']} obj {summ['exact_obj']:.0f} (limit {t}s)", flush=True)
            return {c for c in (canon(C, n) for C in cl) if c}
        except exact_partial.NotCertified as e:
            REPLICATES.append(dict(tag=tag, status=str(e), tl=t, certified=False)); print(f'{tag}: not certified at {t}s, retrying', flush=True); t *= 4
    return None


def compatible(a, b):
    return not (a & b) or a <= b or b <= a


def backbone(stab, theta):
    kept = []
    for C, s in sorted(stab.items(), key=lambda x: (-x[1], -len(x[0]))):
        if s < theta: break
        if all(compatible(C, D) for D in kept): kept.append(C)
    return kept


def parsimony(wit, rows, fams, out, tag, K=200):
    n = len(wit); f = os.path.join(out, f'aug_{tag}.num')
    with open(f, 'w') as h:
        h.write(' '.join(wit) + '\n')
        for r in rows: h.write(' '.join(r) + '\n')
        for C in fams:
            line = ' '.join('1' if i in C else '0' for i in range(n)) + '\n'
            for _ in range(K): h.write(line)
    nwk = os.path.join(out, f'pars_{tag}.nwk')
    subprocess.run(['Rscript', os.path.join(HERE, 'pars.R'), f, nwk], env=RENV, capture_output=True, check=True)
    adj = parse_newick(open(nwk).read(), wit)
    return tree_clusters(adj, n, 0)


def heinrichi_score(wit, fams, truth, out, tag):
    from score import score
    f = os.path.join(out, f'tree_{tag}.txt')
    write_graph(tree_from_clusters(list(fams), len(wit)), wit, f, directed=False)
    return round(score(truth, f)[0], 2)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('num'); ap.add_argument('--out', required=True)
    ap.add_argument('--B', type=int, default=50); ap.add_argument('--thetas', default='0.5,0.7,0.9')
    ap.add_argument('--truth'); ap.add_argument('--tl', type=int, default=300); ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    wit, rows = load_num(a.num); n = len(wit); random.seed(a.seed); t0 = time.time()
    full = mwis_families(wit, rows, a.tl)
    if full is None: sys.exit('full-data optimum not certified')
    counts = {}; ncert = 0
    for b in range(a.B):
        rb = [rows[random.randrange(len(rows))] for _ in rows]
        fams = mwis_families(wit, rb, a.tl, tag=f'rep{b}')
        if fams is None: continue
        ncert += 1
        for C in fams: counts[C] = counts.get(C, 0) + 1
        if (b + 1) % 10 == 0: print(f'bootstrap {b+1}/{a.B} ({time.time()-t0:.0f}s)', flush=True)
    stab = {C: c / ncert for C, c in counts.items()}   # denominator: certified replicates only
    res = dict(full_mwis_families=len(full), B=a.B, certified_replicates=ncert, thetas={})
    score = (lambda fams, tag: heinrichi_score(wit, fams, a.truth, a.out, tag)) if a.truth else (lambda fams, tag: None)
    res['mwis_full_score'] = score(full, 'mwis_full')
    P = parsimony(wit, rows, [], a.out, 'plain'); res['parsimony_plain'] = dict(families=len(P), score=score(P, 'pars_plain'))
    for th in [float(x) for x in a.thetas.split(',')]:
        bb = backbone(stab, th)
        T = parsimony(wit, rows, bb, a.out, f't{th}')
        lost = [C for C in bb if C not in T]
        res['thetas'][str(th)] = dict(backbone=len(bb), backbone_score=score(bb, f'bb{th}'), families=len(T),
                                      score=score(T, f'hyb{th}'), backbone_broken=len(lost))
        print(th, res['thetas'][str(th)], flush=True)
    res['stability'] = sorted([[sorted(wit[i] for i in C), round(s, 2), C in full] for C, s in stab.items() if s >= 0.1],
                              key=lambda x: -x[1])
    res['full_unstable'] = [[sorted(wit[i] for i in C), round(stab.get(C, 0), 2)] for C in full if stab.get(C, 0) < 0.5]
    res['replicates'] = REPLICATES
    json.dump(res, open(os.path.join(a.out, 'protocol.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('stability', 'full_unstable')}))
