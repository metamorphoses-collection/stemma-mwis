"""Q6: classical baselines on the three unweighted graphs run on Ruby, at the hardware's number of copy samples.
(a) random-order greedy maximal independent set: P(result is a maximum independent set), exact over all n! orders;
(b) uniform sampling over all independent sets: share that are maximum, exact;
plus uniform over maximal sets and uniform over all 2^n configurations, and a seeded Monte Carlo draw of exactly the
hardware's number of copy samples for (a) and (b). Output: checks_20261123/baselines.json"""
import itertools as it, math, random
from collections import Counter
from check_common import GRAPHS, LAW_E, LAW_R, conflict_fn, independent_sets, maximum_sets, maximal_sets, lab, dump

HW_N = {'priscus': 31500, 'lai': 8000, 'lawrence': 7000}


def greedy_dist(g):
    conf = conflict_fn(g); F = g['fam']; d = Counter()
    for perm in it.permutations(F):
        s = []
        for v in perm:
            if all(not conf(v, u) for u in s): s.append(v)
        d[frozenset(s)] += 1
    tot = math.factorial(len(F))
    return {k: v / tot for k, v in d.items()}


def greedy_once(g, rng):
    conf = conflict_fn(g); F = list(g['fam']); rng.shuffle(F); s = []
    for v in F:
        if all(not conf(v, u) for u in s): s.append(v)
    return frozenset(s)


out = {}
rng = random.Random(20261123)
for name, g in GRAPHS.items():
    ind = independent_sets(g); mis = maximum_sets(g); mxl = maximal_sets(g); n = len(g['fam'])
    gd = greedy_dist(g)
    p_greedy = sum(gd.get(s, 0) for s in mis)
    p_unif = len(mis) / len(ind)
    N = HW_N[name]
    mc_g = sum(greedy_once(g, rng) in mis for _ in range(N)) / N
    mc_u = sum(rng.choice(ind) in mis for _ in range(N)) / N
    se = lambda p: math.sqrt(p * (1 - p) / N)
    rec = dict(
        vertices=n, edges=len(g['conflicts']), edge_list=g['conflicts'],
        n_independent_sets=len(ind), mis_size=len(mis[0]), n_maximum_sets=len(mis),
        maximum_sets=[lab(g, s) for s in mis],
        maximal_sets=[lab(g, s) for s in mxl], maximal_sizes=sorted(len(s) for s in mxl),
        greedy_P_maximum=p_greedy, greedy_dist={lab(g, s): v for s, v in sorted(gd.items(), key=lambda kv: -kv[1])},
        uniform_IS_P_maximum=p_unif,
        uniform_maximal_P_maximum=sum(s in mis for s in mxl) / len(mxl),
        random_config_P_maximum=len(mis) / 2 ** n,
        equal_shots=dict(N_copy_samples=N, greedy_MC=mc_g, greedy_binom_se=se(p_greedy),
                         uniform_IS_MC=mc_u, uniform_IS_binom_se=se(p_unif)),
        hardware_reported=g['hw_reported'], noisy_emulation_reported=g['emu_noisy_reported'])
    if name == 'lawrence':
        e, r = gd.get(LAW_E, 0), gd.get(LAW_R, 0)
        rec['greedy_E_over_E_plus_R'] = e / (e + r)
        rec['uniform_IS_E_over_E_plus_R'] = 0.5
    out[name] = rec
    print(f"{name}: n={n} m={len(g['conflicts'])} |IS|={len(ind)} MIS size {len(mis[0])} x{len(mis)} | maximal sizes {rec['maximal_sizes']}"
          f" | greedy P(max)={p_greedy:.3f} (MC {mc_g:.3f}, se {se(p_greedy):.4f}) | uniform-IS P(max)={p_unif:.3f} (MC {mc_u:.3f})"
          f" | uniform-maximal {rec['uniform_maximal_P_maximum']:.3f} | random config {rec['random_config_P_maximum']:.4f}"
          f" | hardware {g['hw_reported']} noisy emu {g['emu_noisy_reported']}")
    if name == 'lawrence':
        print('  greedy E/(E+R) =', round(rec['greedy_E_over_E_plus_R'], 4), '; uniform-IS E/(E+R) = 0.5')
print('->', dump(out, 'baselines.json'))
