# -*- coding: utf-8 -*-
"""Fig4/Fig5 v2 — journal-grade redesign.
Fig4 (a) paired dumbbell with KINS-round grouping; (b) individual retentions
+ mean±SD overlay. Fig5 heatmap with Δ±u per cell, significance outlines,
original/replicate divider. Style: Ku journal rules (figsize 20x10, dpi200,
panel labels 28pt y=1.18 ha=left, axis 16 / tick 15, warm palette, no
suptitle) + dataviz craft (recessive axes, thin marks, text in ink colors)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import csv

LABEL, AXIS, TICK = 28, 16, 15
COL = {'Cs-137': '#C94F4A', 'Cs-134': '#E8943A', 'Am-241': '#4AACB0', 'Co-60': '#5B8DB8'}
INK, INK2 = '#2B2B2B', '#6B6560'
plt.rcParams.update({'font.size': TICK, 'axes.labelsize': AXIS,
                     'xtick.labelsize': TICK, 'ytick.labelsize': TICK,
                     'axes.edgecolor': '#B8B2AA', 'axes.linewidth': 1.0,
                     'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.labelcolor': INK, 'text.color': INK})

rows = list(csv.DictReader(open('/home/claude/pectin_rev/out/dataset_retention.csv')))
for r in rows:
    for k in ('A_before','u_b','A_after','u_a','R_ret','uR','Delta'):
        r[k] = float(r[k])
    r['sample'] = int(r['sample'])
def key(r): return r['sample'] if r['set']=='orig' else 5 + r['sample']

# ══════════════ Fig 4 ══════════════
fig, axes = plt.subplots(1, 2, figsize=(20, 10), gridspec_kw={'width_ratios':[1.15,1]})
plt.subplots_adjust(wspace=0.30, left=0.065, right=0.98, top=0.84, bottom=0.11)

# ── (a) paired dumbbell, Cs-137 ──
ax = axes[0]
cs = sorted([r for r in rows if r['nuclide']=='Cs-137'], key=key)
x = np.arange(1, 9)
c = COL['Cs-137']
# KINS round grouping bands
groups = [(0.55, 2.45, '2024 round\n45.1 Bq kg$^{-1}$'),
          (2.55, 4.45, '2022 round\n24.6 Bq kg$^{-1}$'),
          (4.55, 5.45, '2019 round\n24.8 Bq kg$^{-1}$'),
          (5.55, 8.45, '2022 round, replicates\n24.6 Bq kg$^{-1}$')]
for i,(x0,x1,lab) in enumerate(groups):
    if i % 2 == 0:
        ax.axvspan(x0, x1, color='#F4F1EC', zorder=0)
    ax.text((x0+x1)/2, 51.0, lab, ha='center', va='top', fontsize=12.5,
            color=INK2, linespacing=1.25)
ax.set_axisbelow(True)
ax.grid(axis='y', color='#E4E0D9', lw=0.9)
dx = 0.16
for i, r in enumerate(cs):
    xi = x[i]
    # connector
    ax.plot([xi-dx, xi+dx], [r['A_before'], r['A_after']], color=c, lw=1.4,
            alpha=0.45, zorder=2)
    ax.errorbar(xi-dx, r['A_before'], yerr=r['u_b'], fmt='o', ms=10.5,
                mfc='white', mec=c, mew=2.0, ecolor=c, elinewidth=1.4,
                capsize=3.5, zorder=3)
    ax.errorbar(xi+dx, r['A_after'], yerr=r['u_a'], fmt='o', ms=10.5,
                color=c, mec=c, ecolor=c, elinewidth=1.4, capsize=3.5, zorder=3)
ax.set_xticks(x); ax.set_xticklabels([f'S{i}' for i in x])
ax.set_xlim(0.4, 8.6); ax.set_ylim(20, 52)
ax.set_xlabel('Sample')
ax.set_ylabel(r'$^{137}$Cs activity concentration (Bq kg$^{-1}$)')
for s in ('top','right'):
    ax.spines[s].set_visible(False)
handles = [Line2D([0],[0], marker='o', ls='', ms=10.5, mfc='white', mec=c, mew=2.0,
                  label='Before solidification'),
           Line2D([0],[0], marker='o', ls='', ms=10.5, color=c,
                  label='After solidification')]
ax.legend(handles=handles, loc='lower left', fontsize=TICK-1, frameon=False,
          handletextpad=0.4, borderaxespad=0.6)
ax.text(0.0, 1.18, '(a)', transform=ax.transAxes, fontsize=LABEL,
        fontweight='bold', ha='left', va='top')
ax.tick_params(direction='out', length=4.5)

# ── (b) individual retentions + mean±SD ──
ax = axes[1]
order = ['Cs-137', 'Cs-134', 'Am-241', 'Co-60']
ax.axhspan(80, 120, color='#F4F1EC', zorder=0)
for yv in (80, 120):
    ax.axhline(yv, color='#D6D0C7', lw=1.0, zorder=1)
ax.axhline(100, ls=(0,(5,3)), color='#8A847C', lw=1.3, zorder=1)
ax.set_axisbelow(True)
rng = np.random.default_rng(7)
for i, nuc in enumerate(order):
    v = np.array(sorted(r['R_ret'] for r in rows if r['nuclide']==nuc))
    n = len(v)
    jit = (np.arange(n) - (n-1)/2) * (0.30/max(n-1,1))
    ax.scatter(i-0.16+jit, v, s=62, color=COL[nuc], alpha=0.50, lw=0, zorder=3)
    m, sd = v.mean(), v.std(ddof=1)
    ax.errorbar(i+0.17, m, yerr=sd, fmt='D', ms=12, color=COL[nuc], mec='white',
                mew=1.4, ecolor=COL[nuc], elinewidth=2.4, capsize=6,
                capthick=2.4, zorder=4)
    ax.annotate(f'{m:.1f} ± {sd:.1f}', (i+0.17, m), textcoords='offset points',
                xytext=(13, 5), fontsize=13.5, ha='left', va='bottom', color=INK)
ax.text(-0.55, 119.0, '80–120% screening range', fontsize=12, color='#8A847C',
        ha='left', va='top')
ax.set_xticks(range(4))
ax.set_xticklabels([r'$^{137}$Cs' '\n' r'$n=8$', r'$^{134}$Cs' '\n' r'$n=7$',
                    r'$^{241}$Am' '\n' r'$n=8$', r'$^{60}$Co' '\n' r'$n=3$'])
ax.set_xlim(-0.62, 3.62); ax.set_ylim(75, 125)
ax.set_ylabel('Paired retention (%)')
for s in ('top','right'):
    ax.spines[s].set_visible(False)
handles = [Line2D([0],[0], marker='o', ls='', ms=9, color='#9A9C9E', alpha=0.6,
                  label='Individual pairs'),
           Line2D([0],[0], marker='D', ls='', ms=11, color='#6B6560', mec='white',
                  label='Mean ± 1 SD')]
ax.legend(handles=handles, loc='lower right', fontsize=TICK-1, frameon=False,
          handletextpad=0.4, borderaxespad=0.6)
ax.text(0.0, 1.18, '(b)', transform=ax.transAxes, fontsize=LABEL,
        fontweight='bold', ha='left', va='top')
ax.tick_params(direction='out', length=4.5)
fig.savefig('/home/claude/pectin_rev/manuscript/figure/Fig4_revised.png', dpi=200)
plt.close(fig)

# ══════════════ Fig 5 ══════════════
fig, ax = plt.subplots(figsize=(20, 10))
plt.subplots_adjust(left=0.085, right=0.90, top=0.86, bottom=0.10)
grid = np.full((4, 8), np.nan); unc = np.full((4, 8), np.nan); sig = np.zeros((4,8), bool)
for r in rows:
    i, j = order.index(r['nuclide']), key(r)-1
    grid[i, j] = r['Delta']; unc[i, j] = r['uR']
    sig[i, j] = abs(r['Delta']) > 2*r['uR']
cmap = LinearSegmentedColormap.from_list('warmdiv',
        ['#B6423D', '#D98D8A', '#F6E9E8', '#FFFFFF', '#E4F0F1', '#8CC5C8', '#3E8F93'])
im = ax.imshow(grid, cmap=cmap, vmin=-16, vmax=16, aspect='auto', zorder=1)
for i in range(4):
    for j in range(8):
        if np.isnan(grid[i, j]):
            ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, facecolor='#F7F5F1',
                                       zorder=2))
            ax.text(j, i, 'n.m.', ha='center', va='center', fontsize=14,
                    style='italic', color='#B0AAA1', zorder=3)
        else:
            dark = abs(grid[i, j]) >= 10
            c1 = 'white' if dark else INK
            c2 = '#F3E6E5' if dark else INK2
            ax.text(j, i-0.10, f'{grid[i,j]:+.1f}', ha='center', va='center',
                    fontsize=19, color=c1, zorder=3)
            ax.text(j, i+0.24, f'± {unc[i,j]:.1f}', ha='center', va='center',
                    fontsize=12.5, color=c2, zorder=3)
            if sig[i, j]:
                ax.add_patch(plt.Rectangle((j-0.47, i-0.47), 0.94, 0.94,
                             fill=False, edgecolor='#3F3B37', lw=2.0, zorder=4))
# white cell separators
ax.set_xticks(np.arange(-0.5, 8, 1), minor=True)
ax.set_yticks(np.arange(-0.5, 4, 1), minor=True)
ax.grid(which='minor', color='white', lw=3.5)
ax.tick_params(which='minor', length=0)
# original / replicate divider
ax.axvline(4.5, color='white', lw=9, zorder=5)
ax.annotate('Original preparations (S1–S5)', xy=(2.0, -0.62), ha='center',
            va='bottom', fontsize=14.5, color=INK2, annotation_clip=False)
ax.annotate('Replicates, 2022 round (S6–S8)', xy=(6.5, -0.62), ha='center',
            va='bottom', fontsize=14.5, color=INK2, annotation_clip=False)
ax.set_xticks(range(8)); ax.set_xticklabels([f'S{i}' for i in range(1, 9)])
ax.set_yticks(range(4))
ax.set_yticklabels([r'$^{137}$Cs', r'$^{134}$Cs', r'$^{241}$Am', r'$^{60}$Co'],
                   fontsize=17)
ax.set_xlabel('Sample')
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.025)
cb.set_label(r'Paired change $\Delta$ (%)', fontsize=AXIS)
cb.set_ticks([-16, -8, 0, 8, 16])
cb.ax.tick_params(labelsize=TICK, length=0)
cb.outline.set_edgecolor('#B8B2AA'); cb.outline.set_linewidth(0.8)
# significance key
ax.add_patch(plt.Rectangle((7.62, 3.72), 0.22, 0.20, fill=False,
             edgecolor='#5A5550', lw=2.0, clip_on=False, zorder=6))
ax.text(7.92, 3.82, r'$|\Delta| > 2u$', fontsize=13.5, color=INK2, va='center',
        ha='left', clip_on=False, zorder=6)
fig.savefig('/home/claude/pectin_rev/manuscript/figure/Fig5_revised.png', dpi=200)
plt.close(fig)
print('v2 figures written')
