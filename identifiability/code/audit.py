r"""Final audit of the technical report: one entry per result, keyed by its LaTeX label, in report order.

For each result the driver runs the machine check(s) that test it (the existing scripts, unchanged, as subprocesses;
the new checks of audit_new_checks.py in process) and records: label, the report's number (from main.aux when found),
what is tested (one sentence), instances, outcome (pass/fail and the key numbers read from the check's own output)
and the script. The results moved from the mathematical paper (prop:gap, prop:modelbounds, prop:quartet, prop:drift and
the corollaries) are report entries of Section 8, in report order; their former check numbers stay in the JSON.

Outputs: data/final_audit_YYYYMMDD.json and data/final_audit_table.tex (a longtable whose rows are keyed by
\ref{label}, so the report's numbering is used). Exits non-zero if any check fails or any script exits non-zero.
Run from the paper directory (identifiability/):  python3 code/audit.py
Inputs outside this directory: the Lai units (in the MC repository) and the Heinrichi CASC data (~/casc_heinrichi,
not redistributable) for the data part of the alpha check."""
import sys, os, json, re, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); os.chdir(ROOT)
sys.path.insert(0, HERE)
import audit_new_checks as NC

DATE = os.environ.get('AUDIT_DATE', datetime.date.today().strftime('%Y%m%d'))
AUX = os.path.expanduser(os.environ.get('REPORT_AUX', 'main.aux'))
T0 = time.time()

# ------------------------------------------------------------------ existing scripts (the machinery, unchanged)
SCRIPTS = [  # tag, script, args
    ('C1-C14', 'verify_claims', ['data']),
    ('C15', 'kappa_quartet', ['data/kappa_quartet.json']),
    ('C16', 'mechanisms_check', ['data/mechanisms_check.json']),
    ('C17', 'drift_quartet', ['data/drift_quartet.json']),
    ('C27', 'multistate_check', ['data/multistate_check.json']),
    ('C28', 'tworeading_check', ['data/tworeading_check.json']),
    ('C29', 'tops_check', ['data/tops_check.json']),
    ('C18', 'gap_check', ['data/gap_check.json']),
    ('C19', 'fragmentary_check', ['data/fragmentary_check.json']),
    ('C20', 'alpha_estimate', ['data/alpha_estimate.json']),
    ('C21', 'bedier_birthdeath_check', ['data/bedier_birthdeath.json', 'data/bedier_birthdeath_spencer.json',
                                        'data/bedier_birthdeath_spencer_L300.json']),
    ('C22', 'gamma_confidence', ['data/gamma_confidence.json']),
    ('C23', 'polytomy_null', ['data/polytomy_null.json']),
    ('C24', 'bedier_calibrated_check', ['data/bedier_calibrated.json', 'data/bedier_calibrated_spencer.json']),
    ('C25', 'cutting_consistency', ['data/cutting_consistency.json']),
    ('C26', 'cutting_consistency_sim', ['data/cutting_consistency_sim.json']),
    ('rival', 'rival_test', ['data/rival_test_lai.json']),          # computed on the Lai, not a check
    ('rival_unequal', 'rival_unequal', ['data/rival_unequal_lai.json']),   # computed on the Lai, not a check
    ('regime', 'coding_regime', ['data/coding_regime.json']),        # exact boundaries, reported in the paper
]
RUN = {}
for tag, script, args in SCRIPTS:
    t = time.time(); log = f'data/check_{tag}.log'
    with open(log, 'w') as fh:
        rc = subprocess.run([sys.executable, f'code/{script}.py', *args], stdout=fh, stderr=subprocess.STDOUT).returncode
    RUN[tag] = dict(script=script, exit=rc, seconds=round(time.time() - t, 1), log=log)
    print(f'{tag:7s} {script:26s} exit {rc} ({RUN[tag]["seconds"]} s)', flush=True)

NEW = {}
for name in NC.CHECKS:
    NEW[name] = NC.run(name)
json.dump(NEW, open('data/audit_new_checks.json', 'w'), indent=1, default=float)


def J(path): return json.load(open(path))
def LOG(tag): return open(RUN[tag]['log']).read()
def num(x): return f'{x:,}'.replace(',', '{,}')


V = J('data/verify_claims.json')


def from_log(tag):
    t = LOG(tag); p, f = t.count('\nPASS ') + t.startswith('PASS '), t.count('\nFAIL ') + t.startswith('FAIL ')
    return p, f


# ------------------------------------------------------------------ extractors: (instances, result, passed)
def x_C1():
    c = V['C1_universality']; return f'{c["graphs"]} random graphs, $2\\le|V|\\le8$', \
        'all realised ($4|E|+2|V|+2$ witnesses, distinct non-trivial splits, conflict graph equal to $H$)' if c['all_realised'] else 'FAILED', c['all_realised']


def x_C4():
    c = V['C4_maximality']; ok = c['maximal'] and c['max_size_ok']
    return f'{c["trees"]} random binary trees, $4\\le n\\le9$, all splits scanned per tree; sampled compatible sets', \
        'size $n-3$, maximality and the bound $n-3$ hold in all', ok


def x_C2():
    c = V['C2_pairs_triples']; ok = all(c[k]['pairs_completable'] for k in ('4', '5', '6')) and c['4']['obstructing_triples'] == 0 \
        and c['5']['obstructing_triples'] == 200 and c['example_pqr'] == dict(pairwise=True, completable=False)
    return 'every partial split on $n=4,5,6$ witnesses (domains of size $4$ to $n$)', \
        f'all overlap-compatible pairs completable; {c["4"]["obstructing_triples"]} obstructing triples on $4$ witnesses, ' \
        f'{c["5"]["obstructing_triples"]} on $5$; the triple $p,q,r$ pairwise compatible and not completable', ok


def x_C11():
    v = V['C11_proper_partial_triples_n5']; return 'proper partial splits on five witnesses, all triples', \
        f'{v} obstructing triples of quartets (expected $140$)', v == 140


def x_C9():
    c = V['C9_cutting_planes']; return f'{c["instances"]} random weighted families of partial splits, $n=5,6$', \
        f'cut loop reaches the brute-force optimum in all ({c["cuts_total"]} cuts in total)', c['agree']


def x_C3():
    c = V['C3_displacement_matching']; ok = all(v for k, v in c.items() if k != 'trials')
    return f'{num(c["trials"])} sampled pairs $(T,S)$, $4\\le n\\le9$', '$A\\neq\\emptyset$, $|B|\\le|A|$ and the matching exist in all', ok


def x_C10(kind):
    def f():
        c = V['C10_certificates']; return f'{c["instances_" + kind]} random binary trees, $n=5..7$, weights meeting the {kind} hypothesis; ' \
            'all compatible sets enumerated', 'the true set is the unique optimum in all', c['unique_true_optimum']
    return f


def x_C5():
    c = V['C5_rooted']; ok = all(v for k, v in c.items() if k != 'trees')
    return f'{c["trees"]} random rooted trees, $3\\le n\\le8$ ($\\X$ adjoined)', 'counting and maximality facts and the laminar bound hold in all', ok


def x_C8():
    c = V['C8_gamma']; return f'{c["trials"]} instances, $n=6$, all compatible sets enumerated, no solver', \
        f'bisection equals the ratio definition (tolerance $10^{{-9}}$, largest difference ${NC.tex_e(c["max_abs_diff"])}$)', c['agree']


def x_C12_size():
    c = V['C12_structure']; vc = [c['vertex_counts'][str(n)] for n in range(4, 8)]
    ok = vc == [2 ** (n - 1) - n - 1 for n in range(4, 8)]
    return 'full conflict graphs on $n=4..7$ witnesses', f'vertex counts ${",".join(map(str, vc))}$, equal to $2^{{n-1}}-n-1$', ok


def x_C12_five():
    c = V['C12_structure']; ok = c['n5_iso_line_graph_K5'] and c['n5_iso_complement_petersen'] and not c['n5_chordal'] \
        and c['n5_induced_C4'] and c['n5_induced_C5']
    return 'the full conflict graph on five witnesses', 'isomorphic to $L(K_5)$ and to the complement of the Petersen graph; ' \
        'induced $C_4$ and $C_5$, not chordal', ok


def x_C13():
    c = V['C13_crossing_set']; ok = all(v for k, v in c.items() if k not in ('trees', 'false_splits'))
    return f'{c["trees"]} random binary trees, $n=5..9$, all {num(c["false_splits"])} false splits', \
        'crossing set equals $T[A]\\cap T[B]$, non-empty and connected, in all', ok


def x_C14():
    c = V['C14_dynamic_programme']; g = J('data/gamma_confidence.json')
    return f'{c["instances"]} random instances, $n=5..8$ (and $32$ more, $n=7,10$, in the margin check)', \
        'the dynamic programme equals an exact MWIS (CP-SAT) in all', c['agree'] and g['dp_equals_cpsat']


def x_C19():
    d = J('data/fragmentary_check.json')
    return f'{len(d["instances"])} random rooted trees, $n=5$, random rates, $\\alpha$ and independent illegibility; all innovation ' \
        'outcomes and legibility patterns enumerated', \
        f'proposals are restrictions; both bounds of (b) hold; in the {d["refined_condition_cases"]} instances meeting the sufficient ' \
        f'condition the true set is the unique maximum of $\\tilde m$, and in all {d["scalar_forced_cases"]} with $\\varrho>0$, rerun at ' \
        '$\\delta=0.9\\varrho/(4+2\\mu)$', d['passed']


def x_C22():
    d = J('data/gamma_confidence.json'); rows = d['clean']['rows'] + d['noisy']['rows']
    cov = sum(r['coverage'] * r['reps'] for r in rows); tot = sum(r['reps'] for r in rows)
    ok = d['dp_equals_cpsat'] and all(r['coverage'] >= 0.95 and r['certified'] == r['certified_and_true'] for r in rows)
    return f'two ten-witness regimes with exact $m$, $5$ lengths ($L=1{{,}}600$--$409{{,}}600$), $40$ collations each; the \\emph{{Lai}}', \
        f'$\\gamma_{{\\mathrm{{lo}}}}\\le\\gamma$ in {round(cov)} of {tot} collations (at least $95\\%$ required at each length); every ' \
        f'certified $\\hat S_L$ is the true one; \\emph{{Lai}}: $\\hat S_L$ binary and optimal, $\\gamma_{{\\mathrm{{lo}}}}={d["lai"]["gamma_lo"]:.4f}$, ' \
        'increasing with scaled length', ok


def x_C25():
    d = J('data/cutting_consistency.json'); ok = d['equivalence_pass'] and d['bound_pass'] and d['regression_finding11']
    smin = min(r['slack_min'] for r in d['rows'])
    return f'(a) {d["equivalence_tests"]} random weighted collations, $n=5,6$; (b) the five-witness example; (c) {d["bound_instances"]} random ' \
        'rooted trees, $n=5$, all innovation states and legibility patterns', \
        f'(a) the maximum completable family has the weight of the best displaying tree; (b) uniquely best tree with a split never ' \
        f'observed in full; (c) the inequality holds against every rival tree (smallest slack ${smin:.2f}$)', ok


def x_C26():
    d = J('data/cutting_consistency_sim.json'); A = [r for r in d['rows'] if r['case'] == 'A']; B = [r for r in d['rows'] if r['case'] == 'B']
    return f'six random rooted trees, $n=5$, attraction model with random illegibility, $M$ computed exactly; {d["reps"]} collations per length', \
        f'true tree the unique maximum of $M$: selection frequency from ${min(r["freq_true_selected"][0] for r in A):.2f}$--' \
        f'${max(r["freq_true_selected"][0] for r in A):.2f}$ at $L=25$ to ${min(r["freq_true_selected"][-1] for r in A):.2f}$ at $L=1{{,}}600$; ' \
        f'rival with larger $M$: falls to ${max(r["freq_true_selected"][-1] for r in B):.2f}$', d['passed']


def x_C21():
    p, f = from_log('C21')
    return 'index code on hand-made trees; the saved outputs of the three birth--death settings, all $\\alpha$', \
        f'{p} assertions pass, {f} fail: truth inside the null band, raw optimum over-resolves, rival test and thresholds as stated, ' \
        '$P(\\deg\\ge3)$ ranges (orderings and bands only)', f == 0 and p > 0


def x_C24():
    p, f = from_log('C24')
    return 'the saved outputs, three $\\alpha$ and $300$ traditions in each of the two settings', \
        f'{p} assertions on every number and range stated for the pruning rule pass, {f} fail (reconstruction only, no null rejection rate)', f == 0 and p > 0


def x_C16():
    d = J('data/mechanisms_check.json')
    return f'{len(d["trees"])} random rooted trees, $n=5,6$, random rates and $\\alpha$, all outcomes enumerated', \
        'classification of all outcomes with at most two innovations, and the two-sided bounds on every support, hold in all', d['passed']


def x_C15():
    d = J('data/kappa_quartet.json'); c = d['checks']
    return 'exact supports on the quartet (symbolic); $4\\times10^5$ simulated units at a random parameter point', \
        f'reduce to Proposition~\\ref{{prop:quartet}} at $\\alpha=1$; Monte Carlo largest $|z|={c["monte_carlo_max_z"]}$ (fails at $5$)', d['passed']


def x_C20():
    d = J('data/alpha_estimate.json')
    # Printed figure: the guarded shared-departure ratio of Table 2 of the report (pairs isolated among the
    # legible witnesses, data/remainder_table.json), not an estimate of alpha^2.
    g = (J('data/remainder_table.json').get('Heinrichi') or {}) if os.path.exists('data/remainder_table.json') else {}
    lai = d['lai']
    if g.get('R'):
        data = (f'; data: shared-departure ratio $\\hat C/\\hat R$, Heinrichi ${g["C"]}/{g["R"]}={g["C"]/g["R"]:.2f}$ '
                f'({g["pairs_with_eligible_units"]} pairs isolated among the legible witnesses), \\emph{{Lai}} '
                f'${lai["C"]}/{lai["R"]}={lai["C"]/lai["R"]:.2f}$ (computed)')
    else:
        data = '; Heinrichi data absent'
    return f'60 random rooted trees, $n=5..7$, {d["part1_pairs"]} isolated pairs, all two-innovation outcomes', \
        f'$C^{{(2)}}=\\alpha^2R^{{(2)}}$ exactly{data}', d['part1_pass']


def x_C23():
    d = J('data/polytomy_null.json'); r = [x['remainder_ratio_halving'] for x in d['rows']]
    med = sorted(x['rel_remainder_eps0004'] for x in d['rows'])[len(r) // 2]
    return f'{d["instances"]} random rooted trees with a polytomy ($k=3$ or $4$; {d["root_polytomies"]} at the root), all $3^{{|E|}}$ ' \
        'innovation states at $\\varepsilon=0.004$ and $0.002$', \
        f'the remainder divides by ${min(r):.1f}$--${max(r):.1f}$ when $\\varepsilon$ is halved (required $7$--$9$); median relative ' \
        f'size ${100 * med:.1f}\\%$ (output). Not checked: $\\kappa$, the overwrite case at $\\alpha=1$, the estimators', d['all_pass']


def x_C6():
    c = V['C6_quartet_formulas']; return f'{c["trials"]} random rate vectors', 'the three formulas equal exhaustive enumeration ($10^{-12}$)', c['all_equal']


def x_C7():
    c = V['C7_irreversible_bounds']; return f'{c["trials"]} random rooted trees, $n=6,7$, all mark patterns', 'both bounds hold in all', \
        c['true_lower'] and c['false_upper']


def x_C17():
    d = J('data/drift_quartet.json')
    return 'exact supports (symbolic); $3\\times10^5$ simulated units at a random point', \
        f'supports in $[0,1]$; Monte Carlo largest $|z|={d["monte_carlo_max_z"]}$ (fails at $5$)', d['passed']


def x_C27():
    d = J('data/multistate_check.json'); ks = ",".join(str(k) for k in sorted({r['k'] for r in d['trials']}))
    what = str(len(d['trials'])) + ' random rooted trees ($n=5,6$), $k\\in\\{' + ks + '\\}$ attractive readings, all outcomes enumerated'
    return what, 'two-innovation parts equal to those at $\\alpha=\\sqrt A$; $|m-m^{(\\sqrt A)}|\\le\\Pr(|M|\\ge3)$; dominance', d['ok']


def x_C28():
    d = J('data/tworeading_check.json')
    what = str(len(d['trials'])) + ' random rooted trees ($n=5,6,7$), every outcome with at most two innovations; the quartet supports summed over all outcomes'
    return what, 'only disjoint pairs with both readings attractive propose a false split; $m_{13|24}=\\alpha^2p_1p_3(1-r)^2$, $m_{14|23}=0$, $m_{12|34}=2r+O(\\varepsilon^3)$', d['ok']


def x_C29():
    d = J('data/tops_check.json')
    return f'{d["weightings"]} random weightings on $n={d["n"]}$, every optimum ({d["optima"]}), {d["polytomies"]} vertices with three or more children', \
        'no union of two children has positive weight', d['ok']


def x_C18():
    d = J('data/gap_check.json'); return f'{num(d["instances"])} random instances meeting the hypothesis', \
        'the cut separates true from false in all', d['lemma_passed']


def x_new(name):
    def f():
        r = NEW[name]; return r['instances'], r['result'], r['passed']
    return f


def x_rival_lai():
    d = J('data/rival_test_lai.json')
    pv = lambda x: f'{x:.2g}' if x >= 1e-3 else NC.tex_e(x)
    return '', 'on the \\emph{Lai} (computed, not a check), one-sided $p$ against the strongest rival: ' + \
        ', '.join(f'{r["family"]} ${pv(r["p_one_sided"])}$' for r in d), True


# ------------------------------------------------------------------ the results, in report order
P, L_, Th, C, S = 'Prop.', 'Lemma', 'Thm', 'Cor.', 'Section'
NAME = 'audit\\_new\\_checks'
ENTRIES = [
    ('Preliminaries, complexity and fragmentary witnesses (Sections~\\ref{sec:prelim}--\\ref{sec:partial})', None),
    ('prop:buneman', P, "Buneman's theorem", 'a set of non-trivial splits is a stemma\'s split set iff pairwise compatible, the tree unique up to isomorphism',
     [('new', x_new('buneman'), NAME)]),
    ('lem:facts', L_, 'Two counting facts', 'size $n-3$ and maximality of binary split sets; at most $n-3$ splits in any compatible set',
     [('C4', x_C4, 'verify\\_claims')]),
    ('lem:exclusive', L_, 'Exclusivity', 'splits (and clusters) proposed by one unit are pairwise compatible',
     [('new', x_new('exclusive'), NAME)]),
    ('thm:universal', Th, 'Universality', 'the construction realises a given graph as the conflict graph of distinct non-trivial splits',
     [('C1', x_C1, 'verify\\_claims')]),
    ('cor:nphard', C, 'NP-hardness', 'the construction preserves vertex weights and MWIS values, so it is a reduction',
     [('new', x_new('nphard'), NAME)]),
    ('prop:triple', P, 'Pairs exact, triples not', 'overlap-compatible pairs are completable; triples need not be',
     [('C2', x_C2, 'verify\\_claims'), ('C11', x_C11, 'verify\\_claims')]),
    ('thm:cutting', Th, 'Cutting planes', 'the lazy-cut loop returns a maximum jointly completable family',
     [('C9', x_C9, 'verify\\_claims')]),
    ('Certificates of uniqueness (Section~\\ref{sec:certificates})', None),
    ('lem:hall', L_, 'Displacement matching', 'for compatible $S\\neq\\T$: $A\\neq\\emptyset$, $|B|\\le|A|$ and a matching of $B$ into $A$ by conflicts',
     [('C3', x_C3, 'verify\\_claims')]),
    ('thm:local', Th, 'Local separation', 'under the local hypothesis the true set is the unique maximum-weight compatible set',
     [('C10', x_C10('local'), 'verify\\_claims')]),
    ('thm:global', Th, 'Global separation', 'under the global hypothesis the true set is the unique maximum-weight compatible set',
     [('C10', x_C10('global'), 'verify\\_claims')]),
    ('prop:threshold', P, 'Non-binary threshold', 'with false weights $\\le t<$ true weights, the optimum over splits of weight $>t$ is $\\T$',
     [('new', x_new('threshold'), NAME)]),
    ('prop:rooted', P, 'Rooted analogue', 'the counting and maximality facts for clusters',
     [('C5', x_C5, 'verify\\_claims')]),
    ('Consistency and sample size (Section~\\ref{sec:consistency})', None),
    ('thm:pairwise', Th, 'Pairwise bound', 'the bound is at least the error probability of the estimator',
     [('new', x_new('pairwise'), NAME)]),
    ('lem:gamma', L_, 'Computing $\\gamma$', 'bisection on $g$ gives the margin of the ratio definition',
     [('C8', x_C8, 'verify\\_claims')]),
    ('thm:consistency', Th, 'Consistency', '(a) the Bernstein and Hoeffding bounds dominate the error; (b) $\\Pr(\\hat S_L=\\T)\\to0$ when some $m(S)>m(\\T)$',
     [('new', x_new('consistency'), NAME)]),
    ('cor:samplesize', C, 'Sample size', 'the length given by the formula makes the error at most $\\varepsilon$',
     [('new', x_new('samplesize'), NAME)]),
    ('prop:blocks', P, 'Dependent blocks', 'with $J$ i.i.d.\\ blocks and any dependence inside them, $2(n-3+K)e^{-J\\gamma^2/2}$ dominates the error',
     [('new', x_new('blocks'), NAME)]),
    ('thm:gammaconf', Th, 'Confidence bound for $\\gamma$', 'Clopper--Pearson $\\gamma_{\\mathrm{lo}}$ covers $\\gamma$ with probability $\\ge1-\\delta$; no false certificate',
     [('C22', x_C22, 'gamma\\_confidence')]),
    ('prop:gammablocks', P, 'Certifying $\\gamma$ with dependent units', 'with blocks and Hoeffding intervals, $\\gamma_{\\mathrm{lo}}\\le\\gamma$ for all binary stemmata simultaneously with probability $\\ge1-\\delta$',
     [('new', x_new('gammablocks'), NAME)]),
    ('thm:fragmentary', Th, 'Fragmentary witnesses', 'the two bounds of (b) and the sufficient condition of (c), under independent illegibility',
     [('C19', x_C19, 'fragmentary\\_check')]),
    ('prop:textillegible', P, 'Illegibility that depends on the text', 'the bounds of Theorem~\\ref{thm:fragmentary}(b) with $\\delta$, and the condition of (c), when illegibility depends on the readings',
     [('new', x_new('textillegible'), NAME)]),
    ('thm:cuttingconsistency', Th, 'Cutting-plane consistency', 'the cutting-plane estimator maximises $M$ over trees and is consistent exactly when the true tree maximises $M$',
     [('C25', x_C25, 'cutting\\_consistency'), ('C26', x_C26, 'cutting\\_consistency\\_sim')]),
    ('prop:dimfree', P, 'Tolerance independent of $n$', 'the display probability of a false split by a rival tree is at most $2\\max(\\delta,2\\delta^2)$, whatever $n$',
     [('new', x_new('dimfree'), NAME)]),
    ('Non-binary stemmata and coincident error (Section~\\ref{sec:bedier})', None),
    ('prop:bediereffect', P, 'Bédier effect', '$\\Pr(\\hat S_L=\\T)\\le(1-m_\\tau)^L$ (exact part); over-resolution on birth--death traditions',
     [('new', x_new('bediereffect'), NAME), ('C21', x_C21, 'bedier\\_birthdeath\\_check')]),
    ('cor:tops', C, 'Polytomies that survive', 'at a vertex of a rooted optimum with $k\\ge3$ children, no union of two children is proposed',
     [('C29', x_C29, 'tops\\_check')]),
    ('thm:bedier', Th, 'Bédier threshold', 'the Kullback--Leibler bound dominates the error of the thresholded estimator; the second form dominates the first',
     [('new', x_new('bedier'), NAME)]),
    ('prop:gap', P, 'Largest-gap threshold', 'deterministic step: under the hypothesis the largest gap separates true from false',
     [('C18', x_C18, 'gap\\_check')]),
    ('prop:rival', P, 'Rival test', 'given $w(\\sigma)+w(\\tau)=N$, $w(\\sigma)$ is binomial with parameters $N$ and $m_\\sigma/(m_\\sigma+m_\\tau)$',
     [('new', x_new('rival'), NAME), ('computed', x_rival_lai, 'rival\\_test')]),
    ('thm:mechanisms', Th, 'Mechanisms of false support', 'classification of outcomes with at most two innovations; two-sided bounds on every support',
     [('C16', x_C16, 'mechanisms\\_check')]),
    ('cor:dominance', C, 'Comparison of the codings', 'at second order, the rooted support of each cluster of a split not of the tree is at most the unrooted support of the split',
     [('new', x_new('dominance'), NAME)]),
    ('prop:modelbounds', P, 'Support bounds', 'irreversible model, general $n$: the lower bound on true support and the upper bound on false mass',
     [('C7', x_C7, 'verify\\_claims')]),
    ('cor:modelglobal', C, 'A sufficient condition', 'under the hypothesis, global separation holds for the exact supports',
     [('new', x_new('modelglobal'), NAME)]),
    ('prop:quartet', P, 'Quartet', 'the three formulas for the quartet supports of the irreversible model',
     [('C6', x_C6, 'verify\\_claims')]),
    ('cor:lba', C, 'Textual long-branch zone', 'with $p_2=p_4=0$, $r_1=r_2=r$: consistency iff $(1-r)p_1p_3<r(2-p_1-p_3)$, and $13|24$ selected on the reverse side',
     [('new', x_new('lba'), NAME)]),
    ('thm:coding', Th, 'Coding theorem', 'the exact supports of the attraction model on the quartet, both codings',
     [('C15', x_C15, 'kappa\\_quartet')]),
    ('prop:drift', P, 'Drift towards the root', 'the exact supports on $(((1,3),2),4)$, both codings',
     [('C17', x_C17, 'drift\\_quartet')]),
    ('prop:multistate', P, 'Several attractive readings', 'with $A=\\sum_ja_j^2$, the supports are those of the attraction model at $\\alpha=\\sqrt A$ up to $\\Pr(|M|\\ge3)$',
     [('C27', x_C27, 'multistate\\_check')]),
    ('prop:tworeading', P, 'Two-reading coding', 'only coincidence proposes false splits at second order; the quartet supports and boundary',
     [('C28', x_C28, 'tworeading\\_check')]),
    ('prop:alpha', P, 'Estimating $\\alpha$', '$C^{(2)}=\\alpha^2R^{(2)}$ for isolated pairs',
     [('C20', x_C20, 'alpha\\_estimate')]),
    ('prop:polytomynull', P, 'Coincidence null at a polytomy', 'the leading term of $m(\\tau)$ at a polytomy, with an $O(\\varepsilon^3)$ remainder',
     [('C23', x_C23, 'polytomy\\_null')]),
    ('sec:calibrated', S, 'Calibrated Bédier test, pruning rule', 'every number and range stated for the pruning rule, on the saved outputs',
     [('C24', x_C24, 'bedier\\_calibrated\\_check')]),
    ('The graphs of stemmatics (Section~\\ref{sec:graphs})', None),
    ('prop:size', P, 'Size', 'the number of non-trivial splits is $2^{n-1}-n-1$',
     [('C12', x_C12_size, 'verify\\_claims')]),
    ('prop:crossing', P, 'Crossing set', 'the true splits crossed by a false $\\tau=A|B$ are the internal edges of $T[A]\\cap T[B]$, non-empty and connected',
     [('C13', x_C13, 'verify\\_claims')]),
    ('prop:five', P, 'Five witnesses', 'the conflict graph on five witnesses is $L(K_5)$, neither chordal nor perfect',
     [('C12', x_C12_five, 'verify\\_claims')]),
    ('thm:fpt', Th, 'Dynamic programme', 'the $O(3^n)$ dynamic programme returns the maximum-weight compatible set',
     [('C14', x_C14, 'verify\\_claims')]),
    ('Simulations (Section~\\ref{sec:simulations})', None),
    ('sec:birthdeath', S, 'Birth--death traditions', 'the metrics on hand-made trees, and the orderings and bands stated, on the saved outputs',
     [('C21', x_C21, 'bedier\\_birthdeath\\_check')]),
]
SCRIPT_OF = {'C1': 'C1-C14', 'C2': 'C1-C14', 'C3': 'C1-C14', 'C4': 'C1-C14', 'C5': 'C1-C14', 'C6': 'C1-C14', 'C7': 'C1-C14',
             'C8': 'C1-C14', 'C9': 'C1-C14', 'C10': 'C1-C14', 'C11': 'C1-C14', 'C12': 'C1-C14', 'C13': 'C1-C14', 'C14': 'C1-C14',
             'computed': 'rival'}

numbers = {}
if os.path.exists(AUX):
    for lab, n_ in re.findall(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}', open(AUX).read()): numbers[lab] = n_


def evaluate(key, checks):
    parts = []; ok = True
    for cid, fn, script in checks:
        try:
            inst, res, p = fn()
        except Exception as e:
            inst, res, p = '', f'extractor failed: {e!r}', False
        tag = SCRIPT_OF.get(cid, cid)
        if tag in RUN: p = p and (RUN[tag]['exit'] == 0 or tag == 'C1-C14')   # C1-C14 judged per check, script exit below
        parts.append(dict(former_check=cid, instances=inst, result=res, passed=bool(p),
                          script=script.replace('\\_', '_') + '.py'))
        ok &= bool(p)
    return parts, ok


out = dict(date=DATE, report='report_stemma_compatibility/main.tex', results=[], scripts=RUN, new_checks={})
groups = []
for e in ENTRIES:
    if e[1] is None: groups.append(e[0]); out['results'].append(dict(group=e[0])); continue
    label, kind, name, tested, checks = e
    parts, ok = evaluate(label, checks)
    out['results'].append(dict(label=label, number=numbers.get(label), kind=kind, name=name, tested=tested, checks=parts, passed=ok))
out['new_checks'] = {k: dict(passed=v['passed'], seconds=v['seconds'], numbers=v['numbers']) for k, v in NEW.items()}
fails = [r['label'] for r in out['results'] if 'label' in r and not r['passed']] + \
        [f'script {t} exit {r["exit"]}' for t, r in RUN.items() if r['exit'] != 0]
out['failures'] = fails; out['passed'] = not fails; out['seconds'] = round(time.time() - T0, 1)
json.dump(out, open(f'data/final_audit_{DATE}.json', 'w'), indent=1, default=float)


# ------------------------------------------------------------------ LaTeX
def tt(script):
    return '\\texttt{' + script.replace('_', '\\_\\allowbreak ') + '}'


def row(head, tested, parts):
    segs = []
    for p in parts:
        seg = (p['instances'] + ': ' if p['instances'] else '') + p['result']
        segs.append(seg[0].upper() + seg[1:] if seg else seg)
    status = 'pass' if all(p['passed'] for p in parts) else '\\textbf{FAIL}'
    scripts = ', '.join(dict.fromkeys(tt(p['script'][:-3]) for p in parts))
    return f'{head} & {tested} & ' + '. '.join(segs) + f'. \\emph{{({status})}} & {scripts}\\\\\n'


tex = ['% Generated by code/audit.py on ' + DATE + ' from data/final_audit_' + DATE + '.json. Do not edit by hand.\n',
       '{\\small\\setlength{\\tabcolsep}{4pt}\\renewcommand{\\arraystretch}{1.1}\n',
       '\\begin{longtable}{@{}>{\\raggedright\\arraybackslash}p{0.16\\textwidth}>{\\raggedright\\arraybackslash}p{0.23\\textwidth}'
       '>{\\raggedright\\arraybackslash}p{0.42\\textwidth}>{\\raggedright\\arraybackslash}p{0.13\\textwidth}@{}}\n',
       '\\toprule\nResult & What is tested & Instances and result & Script\\\\\n\\midrule\\endhead\n\\bottomrule\\endfoot\n']
for r in out['results']:
    if 'group' in r:
        tex.append(f'\\addlinespace[3pt]\\multicolumn{{4}}{{@{{}}l}}{{\\textbf{{{r["group"]}}}}}\\\\[1pt]\n'); continue
    tex.append(row(f'{r["kind"]}~\\ref{{{r["label"]}}} {r["name"]}', r['tested'], r['checks']))
tex.append('\\end{longtable}}\n')
open('data/final_audit_table.tex', 'w').write(''.join(tex))

print()
for r in out['results']:
    if 'group' in r: continue
    k = r['label']
    print(f'{"PASS" if r["passed"] else "FAIL"}  {k:24s} {numbers.get(k, ""):6s} ' + ', '.join(p['former_check'] for p in r['checks']))
print(f'\n{"FINAL AUDIT PASS" if not fails else "FINAL AUDIT FAILED: " + str(fails)} ({out["seconds"]} s); '
      f'data/final_audit_{DATE}.json, data/final_audit_table.tex')
sys.exit(1 if fails else 0)
