import re, itertools
L={}
for line in open('heinrichi_prompt.txt'):
    m=re.match(r'^(W\d\d): (.*)$', line.strip())
    if m: L[m.group(1)]=m.group(2).replace(' ','')
def alts(s):  # expand 'a/b' alternatives on single characters, e.g. 2346/8 -> 23466 or 23468? treat x/y as either char
    parts=re.split(r'(\w)/(\w)', s)
    out=['']
    i=0
    toks=re.findall(r'\w/\w|[^/]', s)
    for t in toks:
        if '/' in t: out=[o+c for o in out for c in t.split('/')]
        else: out=[o+t for o in out]
    return out
def has(w,s): return any(a in L[w] for a in alts(s))
W=lambda s:['W%02d'%int(x) for x in re.findall(r'\d+',s)]
claims=[
 ('Sonnet','06 07 28','1312133'),('Sonnet','06 07 28','1131144311'),
 ('Sonnet','23 27 32','2471432346/8'),
 ('Sonnet','13 14 17 18 20 36','86115324/5'),('Sonnet','13 14 17 18 20 36','33c/d?23a?91'),
 ('Sonnet','05 16','225143342c'),('Sonnet','37','249143342c'),('Sonnet','05 37','6?63134411'),
 ('Sonnet','02 30','512441241444333215365356347235'),
 ('Sonnet','03','2321823311'),('Sonnet','12','2323824311'),
 ('Sonnet','04 21 29 31 33 34 35','c22b344'),
 ('Mistral','03 04 05 06 07 11 12 15 18 19 21 22 24 25 26 30 31 32 33 34 37','22b'),
 ('Qwen Max','06 07 28','1111'),
]
for who,ws,s in claims:
    g=W(ws); ok=[w for w in g if has(w,s)]; others=[w for w in L if w not in g and has(w,s)]
    print(f"{who:8s} '{s}' in {len(ok)}/{len(g)} named witnesses; also in {len(others)} others")
