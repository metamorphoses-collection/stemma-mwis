r"""Regime of the factor alpha^2 in the coding theorem (thm:coding): exact consistency boundaries (ab)* of the
unrooted and rooted codings on the quartet, rho = 1, p1 = p3 = sqrt(ab) eps, r = eps^2, from the exact supports summed
over the 3^6 edge outcomes (rooted: best laminar family against the true clusters). Usage: coding_regime.py OUT.json"""
import itertools, math
from collections import defaultdict
L=[1,2,3,4]
def supports(r,p1,p3,al):
    edges=['h12','h34','p1','p2','p3','p4']
    rate={'h12':r,'h34':r,'p1':p1,'p2':0,'p3':p3,'p4':0}
    above={1:['p1','h12'],2:['p2','h12'],3:['p3','h34'],4:['p4','h34']}
    mu=defaultdict(float); mr=defaultdict(float)
    for out in itertools.product([0,1,2],repeat=6):
        pr=1.0
        for e,o in zip(edges,out):
            q=rate[e]
            pr*= (1-q) if o==0 else (q*al if o==1 else q*(1-al))
        if pr==0: continue
        o=dict(zip(edges,out))
        read={}
        for w in L:
            v=0
            for e in above[w]:
                if o[e]==1: v='A';break
                if o[e]==2: v=e;break
            read[w]=v
        cls=defaultdict(set)
        for w,v in read.items(): cls[v].add(w)
        S=set(); C=set()
        for v,c in cls.items():
            if 2<=len(c)<=2:
                s=frozenset(c); S.add(min(s,frozenset(L)-s,key=lambda z: sorted(z)))
            if v!=0 and 2<=len(c)<=3: C.add(frozenset(c))
        for s in S: mu[s]+=pr
        for c in C: mr[c]+=pr
    return mu,mr
clusters=[frozenset(c) for k in (2,3) for c in itertools.combinations(L,k)]
def lam(a,b): return not(a&b) or a<=b or b<=a
fams=[]
for k in range(0,4):
    for F in itertools.combinations(clusters,k):
        if all(lam(a,b) for a,b in itertools.combinations(F,2)): fams.append(F)
T=(frozenset({1,2}),frozenset({3,4}))
def margin(ab,eps,al,coding):
    r=eps**2; p=math.sqrt(ab)*eps
    mu,mr=supports(r,p,p,al)
    if coding=='u':
        return mu[frozenset({1,2})]-max(mu[frozenset({1,3})],mu[frozenset({1,4})])
    vT=sum(mr[c] for c in T)
    best=max(sum(mr[c] for c in F) for F in fams if set(F)!=set(T))
    return vT-best
def boundary(eps,al,coding):
    lo,hi=1e-3,1e6
    for _ in range(80):
        mid=math.sqrt(lo*hi)
        if margin(mid,eps,al,coding)>0: lo=mid
        else: hi=mid
    return lo
import json, sys
rows = []
for al,eps in [(0.05,1e-3),(0.05,1e-2),(0.5,1e-3),(0.5,1e-2),(0.5,0.05),(0.5,0.1),(0.5,0.3),(0.95,0.01),(0.95,0.05),(0.95,0.2),(0.2,0.01),(0.2,0.02)]:
    bu=boundary(eps,al,'u'); br=boundary(eps,al,'r')
    rows.append(dict(alpha=al, eps=eps, unrooted=bu, unrooted_rel=bu/2, rooted=br, rooted_rel=br/(2/al**2))); print(rows[-1])
json.dump(rows, open(sys.argv[1], "w"), indent=1)
