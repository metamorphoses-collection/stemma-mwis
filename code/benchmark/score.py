"""Re-implementation of the CASC 2007 'average sign distance' criterion (skeleton + orientation),
written from the logic of the organisers' tripled.py (read, not executed).
Graph file format: n, then n labels, then an n x n matrix; -1 = no edge, else edge length;
an edge i->j with edge[j][i] == -1 is directed, both != -1 is undirected.
Labels starting with a letter are observed witnesses; others are hidden nodes."""
import sys
from collections import deque


def readgraph(fname):
    with open(fname) as f:
        toks = f.read().split('\n')
    n = int(toks[0])
    lab = [toks[1 + i].strip() for i in range(n)]
    edge = {}
    for i in range(n):
        row = toks[1 + n + i].split()
        for j in range(n):
            edge[lab[i], lab[j]] = float(row[j])
    return n, lab, edge


def shortest(n, lab, edge):
    """Shortest paths ignoring direction (min of the two directions), BFS-relaxation as original."""
    path = {}
    for i in range(n):
        for j in range(n):
            path[lab[i], lab[j]] = -1
        q = deque([(lab[i], 0.0)])
        while q:
            cur, d = q.popleft()
            if path[lab[i], cur] == -1 or path[lab[i], cur] > d:
                path[lab[i], cur] = d
                for j in range(n):
                    a, b = edge[cur, lab[j]], edge[lab[j], cur]
                    if a == -1:
                        m = b
                    elif b == -1:
                        m = a
                    else:
                        m = min(a, b)
                    if m != -1 and lab[j] != cur:
                        q.append((lab[j], d + m))
    return path


def roots(n, lab, edge):
    return [lab[i] for i in range(n)
            if not any(i != j and edge[lab[j], lab[i]] != -1 for j in range(n))]


def depth(n, lab, edge, rts):
    dep = {l: -1 for l in lab}
    q = deque([(r, 0.0) for r in rts])
    while q:
        cur, d = q.popleft()
        if dep[cur] == -1 or dep[cur] > d:
            dep[cur] = d
            for j in range(n):
                m = edge[cur, lab[j]]
                if m != -1 and lab[j] != cur:
                    q.append((lab[j], d + m))
    return dep


def sign(x, y):
    if x == -1 and y == -1:
        return 0
    if x == -1:
        return -1
    if y == -1:
        return 1
    return -1 if x < y else (1 if x > y else 0)


def score(truefile, propfile):
    nt, lt, et = readgraph(truefile)
    npr, lp, ep = readgraph(propfile)
    pt, pp = shortest(nt, lt, et), shortest(npr, lp, ep)
    obs = [l for l in lt if l[0].isalpha()]
    s1 = t1 = 0.0
    for i in obs:
        for j in obs:
            for k in obs:
                if i != j and i != k and j != k:
                    st = sign(pt[i, j], pt[i, k])
                    sp = sign(pp.get((i, j), -1), pp.get((i, k), -1))
                    s1 += 1 - abs(st - sp) / 2.0
                    t1 += 1
    dt = depth(nt, lt, et, roots(nt, lt, et))
    dp = depth(npr, lp, ep, roots(npr, lp, ep))
    s2 = t2 = 0.0
    for i in obs:
        for j in obs:
            if i != j:
                s2 += 1 - abs(sign(dt[i], dt[j]) - sign(dp.get(i, -1), dp.get(j, -1))) / 2.0
                t2 += 1
    return 100 * s1 / t1, 100 * s2 / t2


if __name__ == '__main__':
    a, b = score(sys.argv[1], sys.argv[2])
    print(f'skeleton {a:.1f}%  orientation {b:.1f}%')
