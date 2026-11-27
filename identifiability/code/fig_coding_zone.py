r"""Figure: the textual long-branch zone in the two codings (Theorem coding).
Quartet 12|34, p2 = p4 = 0, halves of rate r, pendant rates p1 = p3 = p. For each r, the exact boundary product p^2
at which the false split 13|24 (unrooted) or a family containing {1,3} (rooted) starts to win is found by bisection
on the exact supports summed over the 3^6 edge outcomes. Usage: fig_coding_zone.py OUT.pdf"""
import sys, itertools, math
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

L = [1, 2, 3, 4]


def supports(r, p, al):
    edges = ['h12', 'h34', 'p1', 'p2', 'p3', 'p4']
    rate = {'h12': r, 'h34': r, 'p1': p, 'p2': 0, 'p3': p, 'p4': 0}
    above = {1: ['p1', 'h12'], 2: ['p2', 'h12'], 3: ['p3', 'h34'], 4: ['p4', 'h34']}
    mu = defaultdict(float); mr = defaultdict(float)
    for out in itertools.product([0, 1, 2], repeat=6):
        pr = 1.0
        for e, o in zip(edges, out):
            q = rate[e]; pr *= (1 - q) if o == 0 else (q * al if o == 1 else q * (1 - al))
        if pr == 0: continue
        o = dict(zip(edges, out)); cls = defaultdict(set)
        for w in L:
            v = 0
            for e in above[w]:
                if o[e] == 1: v = 'A'; break
                if o[e] == 2: v = e; break
            cls[v].add(w)
        S, C = set(), set()
        for v, c in cls.items():
            if len(c) == 2:
                s = frozenset(c); S.add(s if 1 in s else frozenset(L) - s)
            if v != 0 and 2 <= len(c) <= 3: C.add(frozenset(c))
        for s in S: mu[s] += pr
        for c in C: mr[c] += pr
    return mu, mr


CL = [frozenset(c) for k in (2, 3) for c in itertools.combinations(L, k)]
lam = lambda a, b: not (a & b) or a <= b or b <= a
FAMS = [F for k in range(0, 3) for F in itertools.combinations(CL, k) if all(lam(a, b) for a, b in itertools.combinations(F, 2))]
T = {frozenset({1, 2}), frozenset({3, 4})}


def margin(r, p, al, coding):
    mu, mr = supports(r, p, al)
    if coding == 'u':
        return mu[frozenset({1, 2})] - max(mu[frozenset({1, 3})], mu[frozenset({1, 4})])
    return sum(mr[c] for c in T) - max(sum(mr[c] for c in F) for F in FAMS if set(F) != T)


def boundary(r, al, coding):
    lo, hi = 1e-6, 0.999
    if margin(r, hi, al, coding) > 0: return None
    for _ in range(40):
        mid = math.sqrt(lo * hi)
        if margin(r, mid, al, coding) > 0: lo = mid
        else: hi = mid
    return lo


rs = np.geomspace(2e-4, 2e-2, 14)
curves = {'unrooted': [boundary(r, 0.6, 'u') for r in rs]}
ALS = {'Heinrichi': 0.82 ** 0.5, 'Donatus': 0.75 ** 0.5, 'Lai': 0.34 ** 0.5, 'illustration': 0.3}
for nm, al in ALS.items():
    curves[nm] = [boundary(r, al, 'r') for r in rs]

INK, INK2, SURF = '#0b0b0b', '#52514e', '#fcfcfb'
STY = {'Heinrichi': ('#86b6ef', '-'), 'Donatus': ('#5598e7', '-'), 'Lai': ('#1c5cab', '-'), 'illustration': ('#104281', (0, (4, 2)))}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'font.family': 'serif'})
fig, ax = plt.subplots(figsize=(5.2, 4.3))
ax.set_xscale('log'); ax.set_yscale('log')
x = rs
yu = np.array([b * b for b in curves['unrooted']])
ax.plot(x, yu, color=INK2, lw=2, zorder=3)
ax.lines[-1].set_label('counting every reading')
LAB = {'Heinrichi': 'shared departures, Heinrichi (ratio 0.82)', 'Donatus': 'shared departures, Donatus (ratio 0.75)',
       'Lai': 'shared departures, \\emph{Lai} ($\\alpha^2=0.34$)', 'illustration': 'shared departures, $\\alpha^2=0.09$'}
LAB['Lai'] = 'shared departures, Lai (ratio 0.34)'
for nm in ('illustration', 'Lai', 'Donatus', 'Heinrichi'):
    y = np.array([np.nan if b is None else b * b for b in curves[nm]])
    c, ls = STY[nm]
    ax.plot(x, y, color=c, lw=2, ls=ls, zorder=3, label=LAB[nm])
ax.text(x[0] * 1.15, 0.6, 'above a curve: the two corrupted\ncopies are grouped (long-branch zone)', ha='left', va='top', fontsize=8, color=INK2)
h, l = ax.get_legend_handles_labels()
order = [4, 3, 2, 1, 0]
ax.legend([h[i] for i in order], [l[i] for i in order], loc='upper center', bbox_to_anchor=(0.45, -0.17), ncol=2, fontsize=7.5, frameon=False)
ax.set_xlim(rs[0], rs[-1]); ax.set_ylim(3e-4, 1)
ax.set_xlabel('rate of innovation on the internal edge, $r$ (each half)')
ax.set_ylabel('product of the two corruption rates, $p_1p_3$')
ax.grid(True, which='major', color='#e8e8e5', lw=0.6)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(sys.argv[1]); print('saved', sys.argv[1])
import json
json.dump({'r': list(map(float, rs)), 'unrooted': curves['unrooted'], **{f'rooted_{k}': curves[k] for k in ALS}},
          open(sys.argv[1].replace('.pdf', '.json'), 'w'), indent=1)
