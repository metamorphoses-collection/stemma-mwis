import re
L={}
for line in open('heinrichi_prompt.txt'):
    m=re.match(r'^(W\d\d): (.*)$', line.strip())
    if m: L[m.group(1)]=m.group(2).replace(' ','')
N=len(next(iter(L.values())))
def seg(w,a,b): return L[w][a-1:b]          # units are 1-based, inclusive
def lacuna(group,a,b):
    allq=[w for w in group if set(seg(w,a,b))=={'?'}]
    others=[w for w in L if w not in group and set(seg(w,a,b))=={'?'}]
    return f"{len(allq)}/{len(group)} of the group lacunose; outside the group lacunose: {len(others)} {others[:8]}"
def reading(group,a,b,s):
    hit=[w for w in group if seg(w,a,b)==s]; out=[w for w in L if w not in group and seg(w,a,b)==s]
    return f"{len(hit)}/{len(group)} of the group read '{s}'; outside: {len(out)} {out[:8]}"
G=lambda s:['W%02d'%int(x) for x in re.findall(r'\d+',s)]
B=G('08 09 10 13 14 15 17 18 19 20 22 23 26 27 36'); Gb=G('01 04 21 29 31 33 34 35'); H=G('02 03 05 11 12 16 24 25 30 32 37')
claims=[
 ("H branch lacuna 701-711",lambda:lacuna(H,701,711)),
 ("H branch lacuna 813-819",lambda:lacuna(H,813,819)),
 ("F (03 12 11 25 32) lacuna 163-165",lambda:lacuna(G('03 12 11 25 32'),163,165)),
 ("J (02 30 05 16 37) lacuna 665-667",lambda:lacuna(G('02 30 05 16 37'),665,667)),
 ("G branch lacuna 270-272",lambda:lacuna(Gb,270,272)),
 ("G sub (04 21 34 35) lacuna 347-348",lambda:lacuna(G('04 21 34 35'),347,348)),
 ("G sub (04 21 34 35) lacuna 1202-1206",lambda:lacuna(G('04 21 34 35'),1202,1206)),
 ("B branch lacuna 322-325",lambda:lacuna(B,322,325)),
 ("B branch lacuna 1220-1222",lambda:lacuna(B,1220,1222)),
 ("B branch reading '?331' at 541-544",lambda:reading(B,541,544,'?331')),
 ("B minus 17,20 lacuna 845-847",lambda:lacuna([w for w in B if w not in ('W17','W20')],845,847)),
 ("G+B reading '2223' at 121-124",lambda:reading(Gb+B,121,124,'2223')),
 ("F/W24 '353' at 665-667",lambda:reading(G('03 12 11 25 32 24'),665,667,'353')),
]
print('units',N)
for name,f in claims:
    try: print(f"{name:40s} -> {f()}")
    except Exception as e: print(name,'ERR',e)
print('\n-- offset search (does the claim hold if units are shifted by -12..+12?)')
def ok_lac(group,a,b): return all(set(seg(w,a,b))=={'?'} for w in group)
def ok_read(group,a,b,s): return all(seg(w,a,b)==s for w in group)
for name,grp,a,b,s in [("H 701-711",H,701,711,None),("H 813-819",H,813,819,None),("G sub 347-348",G('04 21 34 35'),347,348,None),("G sub 1202-1206",G('04 21 34 35'),1202,1206,None),("B 1220-1222",B,1220,1222,None),("F+W24 '353' 665-667",G('03 12 11 25 32 24'),665,667,'353')]:
    hits=[d for d in range(-12,13) if 1<=a+d and b+d<=N and (ok_lac(grp,a+d,b+d) if s is None else ok_read(grp,a+d,b+d,s))]
    print(f"{name:24s} holds at offsets: {hits}")
