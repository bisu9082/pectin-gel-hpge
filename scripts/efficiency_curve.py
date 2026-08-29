# -*- coding: utf-8 -*-
"""Empirical efficiency calibration from the certified aqueous standard.

Emission rates are those of the KINS-supplied Eckert & Ziegler multi-nuclide
standard (source ID 133399, reference 2025-04-01 12:00), decay-corrected to
each acquisition. Produces the measured full-energy-peak efficiency at each
fill height (2-40 mm) shown in Supplementary Figure S1.

Outputs: raw_data/efficiency_measured.json, figures/FigS1_efficiency.png
"""
import json, datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _paths import SPECTRA, RAW, FIG
from cnf_reader import read_cnf
from spectrum_analysis import peak_net

REF = datetime.datetime(2025, 4, 1, 12, 0)
# nuclide, energy (keV), emission rate at reference date (gamma/s), half-life (d)
GPS = [('Am-241', 59.54, 25.7097, 432.6*365.25), ('Cd-109', 88.03, 80.0656, 461.9),
       ('Co-57', 122.06, 42.4703, 271.8), ('Ce-139', 165.86, 60.3765, 137.64),
       ('Cr-51', 320.08, 467.964, 27.70), ('Sn-113', 391.69, 84.9706, 115.09),
       ('Sr-85', 514.01, 173.035, 64.85), ('Cs-137', 661.66, 53.7755, 30.08*365.25),
       ('Y-88', 898.04, 209.349, 106.63), ('Co-60', 1173.23, 100.377, 5.2711*365.25),
       ('Co-60', 1332.49, 100.389, 5.2711*365.25), ('Y-88', 1836.06, 221.638, 106.63)]
HEIGHTS = [(2, 'PM1EFF255D02'), (5, 'PM1EFF255D05'), (10, 'PM1EFF255D10'),
           (20, 'PM1EFF255D20'), (30, 'PM1EFF255D30'), (40, 'PM1EFF255D40')]
CASCADE = {'Y-88', 'Co-60'}
INK, INK2 = '#2B2B2B', '#6B6560'


def main():
    FIG.mkdir(exist_ok=True)
    out = {}
    for h, stem in HEIGHTS:
        d = read_cnf(SPECTRA/'efficiency_calibration'/f'{stem}.CNF')
        dt = (d['start']-REF).total_seconds()/86400.0
        res = []
        for nuc, E, gps, t12 in GPS:
            r = peak_net(d['counts'], d['ecal'], E)
            decay = 2**(-dt/t12)
            res.append([nuc, E, round(r['net']/d['livetime']/(gps*decay), 5),
                        round(r['u_net']/d['livetime']/(gps*decay), 5)])
        out[h] = dict(days_since_reference=round(dt, 1), efficiency=res)
    json.dump(out, open(RAW/'efficiency_measured.json', 'w'), indent=1)

    plt.rcParams.update({'font.size': 15, 'axes.labelsize': 16,
                         'xtick.labelsize': 15, 'ytick.labelsize': 15,
                         'axes.edgecolor': '#B8B2AA', 'axes.linewidth': 1.0,
                         'xtick.color': INK2, 'ytick.color': INK2,
                         'axes.labelcolor': INK, 'text.color': INK})
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    plt.subplots_adjust(wspace=0.32, left=0.07, right=0.97, top=0.84, bottom=0.12)

    ax = axes[0]
    d10 = out[10]['efficiency']
    E = np.array([r[1] for r in d10]); ep = np.array([r[2] for r in d10])
    ue = np.array([r[3] for r in d10])
    mask = np.array([r[0] not in CASCADE for r in d10])
    fit = mask & (E >= 88)
    c = np.polyfit(np.log(E[fit]), np.log(ep[fit]), 3)
    Eg = np.geomspace(86, 720, 150)
    ax.plot(Eg, np.exp(np.polyval(c, np.log(Eg))), color='#C9C3BA', lw=2.2, zorder=1)
    ax.errorbar(E[mask], ep[mask], yerr=ue[mask], fmt='o', ms=11, color='#C94F4A',
                ecolor='#C94F4A', capsize=4, lw=0, elinewidth=1.5, zorder=3,
                label='Single-line nuclides')
    ax.errorbar(E[~mask], ep[~mask], yerr=ue[~mask], fmt='o', ms=11, mfc='white',
                mec='#5B8DB8', mew=2.2, ecolor='#5B8DB8', capsize=4, lw=0,
                elinewidth=1.5, zorder=3, label='Cascade nuclides (Y-88, Co-60)')
    for nuc, en, e_, u_ in d10:
        ax.annotate(nuc, (en, e_), textcoords='offset points',
                    xytext=(7, 7) if nuc != 'Cd-109' else (-12, -18),
                    fontsize=12, color=INK2)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(45, 2400); ax.set_ylim(0.009, 0.16)
    ax.set_xlabel('Gamma energy (keV)')
    ax.set_ylabel('Full-energy-peak efficiency')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.legend(loc='lower left', fontsize=14, frameon=False)
    ax.text(0.0, 1.18, '(a)', transform=ax.transAxes, fontsize=28,
            fontweight='bold', ha='left', va='top')

    ax = axes[1]
    hs = [h for h, _ in HEIGHTS]
    for idx, colr, lab in [(0, '#4AACB0', '59.5 keV (Am-241)'),
                           (7, '#C94F4A', '661.7 keV (Cs-137)'),
                           (9, '#5B8DB8', '1173.2 keV (Co-60)')]:
        y = [out[h]['efficiency'][idx][2] for h in hs]
        u = [out[h]['efficiency'][idx][3] for h in hs]
        ax.errorbar(hs, y, yerr=u, fmt='o-', ms=10, color=colr, lw=1.6,
                    capsize=4, elinewidth=1.3, label=lab)
    ax.axvline(10, color='#8A847C', ls=(0, (5, 3)), lw=1.3)
    ax.text(10.6, 0.118, 'fill height used\nin this study', fontsize=12.5,
            color='#8A847C', va='top')
    ax.set_xlabel('Fill height (mm)')
    ax.set_ylabel('Full-energy-peak efficiency')
    ax.set_xlim(0, 43); ax.set_ylim(0, 0.125)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.legend(loc='upper right', fontsize=14, frameon=False)
    ax.text(0.0, 1.18, '(b)', transform=ax.transAxes, fontsize=28,
            fontweight='bold', ha='left', va='top')
    fig.savefig(FIG/'FigS1_efficiency.png', dpi=200)
    print(f'Written: {RAW/"efficiency_measured.json"}, {FIG/"FigS1_efficiency.png"}')
    print('Measured efficiency at 10 mm: ' + ', '.join(
        f'{r[0]} {r[1]:.1f} keV = {r[2]:.4f}' for r in d10[:1] + d10[7:8] + d10[9:10]))


if __name__ == '__main__':
    main()
