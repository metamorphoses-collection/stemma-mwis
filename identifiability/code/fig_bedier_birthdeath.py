r"""Figure: the Bédier effect on birth-death traditions (Section 5.2). Share of two-branched roots, i_3 and share of
internal polytomies for the true stemmata and for best-supported reconstructions with decreasing support thresholds,
alpha = 0.6, three settings, against the editors' figures of Camps, Randon-Furling and Godreau.
Reads data/bedier_birthdeath*.json. Usage: fig_bedier_birthdeath.py OUT.pdf"""
import sys, json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D = os.path.join(os.path.dirname(__file__), '..', 'data')
SET = [('equal rates, $L=1{,}000$', 'bedier_birthdeath.json', '#2a78d6'),
       ('gamma rates, $L=1{,}000$', 'bedier_birthdeath_spencer.json', '#eb6834'),
       ('gamma rates, $L=300$', 'bedier_birthdeath_spencer_L300.json', '#1baf7a')]
KEYS = ['true', 'k>=5', 'k>=3', 'k>=2', 'raw']
XL = ['true\nstemma', '$\\geq5$', '$\\geq3$', '$\\geq2$', 'none']
XP = [0, 1.6, 2.6, 3.6, 4.6]
PANELS = [('bifid_root', 'two-branched roots', (0.65, 0.86), None, 'share of stemmata whose\nroot has two branches'),
          ('i3', 'resolved triples, $i_3$', (0.89, 0.96), None, 'share of three-witness subtrees\nresolved into a pair and one'),
          ('p_deg_ge3', 'internal polytomies', None, None, 'share of internal nodes with\nthree or more branches')]
INK, INK2, BAND = '#0b0b0b', '#52514e', '#e8e8e5'
plt.rcParams.update({'font.size': 8.5, 'font.family': 'serif', 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.2), sharey=False)
for ax, (key, title, band, point, ylab) in zip(axes, PANELS):
    ax.set_ylabel(ylab, fontsize=7.5, color=INK2)
    if band: ax.axhspan(*band, color=BAND, zorder=0, lw=0); ax.text(0.8, band[1] - 0.004, "editors'\ninterval", ha='center', va='top', fontsize=7, color=INK2)
    if point is not None: ax.axhline(point, color=INK2, lw=1, ls=(0, (3, 2)), zorder=1); ax.text(0.8, point + 0.008, "editors'", ha='center', va='bottom', fontsize=7, color=INK2)
    ax.axvline(0.8, color=BAND, lw=0.8, zorder=0)
    for j, (name, f, c) in enumerate(SET):
        d = json.load(open(os.path.join(D, f)))['alpha=0.6']
        v = [d[k][key] for k in KEYS]
        ax.plot(XP[1:], v[1:], color=c, lw=2, marker='o', ms=4.5, mec='white', mew=1, label=name, zorder=3)
        ax.plot([XP[0] + (j - 1) * 0.22], [v[0]], color=c, marker='D', ms=4.5, mec='white', mew=1, ls='none', zorder=3)
    ax.set_xticks(XP); ax.set_xticklabels(XL, fontsize=7.5)
    ax.set_title(title, fontsize=9, color=INK)
    ax.set_xlim(-0.4, 4.9)
    ax.grid(True, axis='y', color=BAND, lw=0.6)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
axes[0].set_ylim(0.4, 1.03)
axes[0].axhline(105 / 110, color=INK2, lw=1, ls=(0, (1, 1.5)), zorder=1)
axes[0].text(0.8, 105 / 110 + 0.008, 'Bédier 105/110; MGH 90/94', ha='left', va='bottom', fontsize=7, color=INK2); axes[1].set_ylim(0.6, 1.0); axes[2].set_ylim(0, 0.45)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.02))
fig.supxlabel('true stemma, or best-supported reconstruction admitting only families proposed by at least 5, 3 or 2 units, or by any', fontsize=8, color=INK2, y=0.1)
fig.tight_layout(rect=(0, 0.13, 1, 1))
fig.savefig(sys.argv[1]); print('saved', sys.argv[1])
