# -*- coding: utf-8 -*-
"""Independent re-analysis of all 26 paired spectra.

For each spectrum and each gamma line: net peak area by side-band background
subtraction, counting uncertainty, and Currie critical level / detection limit.
Then paired retention per line (efficiency cancels in the ratio) for the
primary/secondary-line consistency check reported in Supplementary Table S2.

Outputs (raw_data/): peaks_all.csv, recovery_lines.csv
"""
import csv, math, datetime
import numpy as np
from _paths import SPECTRA, RAW
from cnf_reader import read_cnf

LINES = [('Am-241', 59.54, 0.359), ('Cs-134', 604.72, 0.976),
         ('Cs-134', 795.86, 0.855), ('Cs-137', 661.66, 0.851),
         ('Co-60', 1173.23, 0.9985), ('Co-60', 1332.49, 0.9998)]
HALF_LIFE_D = {'Cs-137': 30.08*365.25, 'Cs-134': 2.0652*365.25,
               'Am-241': 432.6*365.25, 'Co-60': 5.2711*365.25}
# sample masses (g) as entered in the Genie-2000 analysis
MASS = {1: (28.2664, 28.3941), 2: (28.2692, 28.7904), 3: (28.2665, 29.1465),
        4: (28.2675, 29.2001), 5: (28.2664, 29.0141), 6: (28.2685, 28.2685),
        7: (28.2683, 28.2683), 8: (28.2665, 28.2665)}


def ch_of_energy(E, cc, nch):
    ch = np.arange(nch)
    en = cc[0] + cc[1]*ch + cc[2]*ch**2 + cc[3]*ch**3
    return float(np.interp(E, en, ch))


def peak_net(counts, cc, E, pk_hw=3.5, bg_gap=4.5, bg_w=4.5):
    """Net area, uncertainty and Currie limits for one line."""
    nch = len(counts)
    c = lambda e: int(round(ch_of_energy(e, cc, nch)))
    p1, p2 = c(E-pk_hw), c(E+pk_hw)
    l1, l2 = c(E-bg_gap-bg_w), c(E-bg_gap)
    r1, r2 = c(E+bg_gap), c(E+bg_gap+bg_w)
    G = counts[p1:p2+1].sum()
    n_pk = p2-p1+1
    BL, BR = counts[l1:l2+1].sum(), counts[r1:r2+1].sum()
    scale = n_pk/((l2-l1+1)+(r2-r1+1))
    B = scale*(BL+BR)
    return dict(gross=float(G), bg=float(B), net=float(G-B),
                u_net=float(math.sqrt(G + B*scale)),
                Lc=float(2.326*math.sqrt(2*B*scale)) if B > 0 else 2.326,
                Ld=float(2.71 + 4.65*math.sqrt(B*scale)) if B > 0 else 2.71)


def main():
    spec, rows = {}, []
    for s in range(1, 9):
        for phase in ('before', 'after'):
            fn = SPECTRA / phase / f'S{s}_{phase}.CNF'
            d = read_cnf(fn)
            spec[(s, phase)] = d
            for nuc, E, br in LINES:
                r = peak_net(d['counts'], d['ecal'], E)
                rows.append(dict(sample=s, phase=phase, nuclide=nuc, E=E,
                                 start=d['start'].isoformat(),
                                 livetime=d['livetime'], **r,
                                 rate=r['net']/d['livetime'],
                                 u_rate=r['u_net']/d['livetime'],
                                 detected=r['net'] > r['Lc']))
    with open(RAW/'peaks_all.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    get = lambda s, ph, nuc, E: next(r for r in rows if (r['sample'], r['phase'],
                                    r['nuclide'], r['E']) == (s, ph, nuc, E))
    rec = []
    for s in range(1, 9):
        mb, ma = MASS[s]
        dt_d = (spec[(s, 'after')]['start'] - spec[(s, 'before')]['start']).total_seconds()/86400
        for nuc, E, br in LINES:
            b, a = get(s, 'before', nuc, E), get(s, 'after', nuc, E)
            if b['net'] <= b['Lc'] or a['net'] <= a['Lc']:
                rec.append(dict(sample=s, nuclide=nuc, E=E, R=None, uR=None,
                                note='below Currie critical level'))
                continue
            decay = math.exp(math.log(2)/HALF_LIFE_D[nuc]*dt_d)
            R = a['rate']/b['rate']*decay*100.0            # retention (mass cancels)
            uR = R*math.sqrt((a['u_rate']/a['rate'])**2 + (b['u_rate']/b['rate'])**2)
            rec.append(dict(sample=s, nuclide=nuc, E=E, R=round(R, 2),
                            uR=round(uR, 2), note=''))
    with open(RAW/'recovery_lines.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['sample', 'nuclide', 'E', 'R', 'uR', 'note'])
        w.writeheader(); w.writerows(rec)

    print('Primary vs secondary line consistency (Supplementary Table S2)')
    for s in range(1, 9):
        for nuc, e1, e2 in [('Cs-134', 604.72, 795.86), ('Co-60', 1173.23, 1332.49)]:
            r1 = next((r for r in rec if r['sample'] == s and r['nuclide'] == nuc
                       and r['E'] == e1 and r['R']), None)
            r2 = next((r for r in rec if r['sample'] == s and r['nuclide'] == nuc
                       and r['E'] == e2 and r['R']), None)
            if r1 and r2:
                z = (r1['R']-r2['R'])/math.sqrt(r1['uR']**2 + r2['uR']**2)
                print(f"  S{s} {nuc:7s} {e1:7.1f}: {r1['R']:6.1f}+-{r1['uR']:4.1f} | "
                      f"{e2:7.1f}: {r2['R']:6.1f}+-{r2['uR']:4.1f}  z={z:+.2f}")
    print(f'\nWritten: {RAW/"peaks_all.csv"}, {RAW/"recovery_lines.csv"}')


if __name__ == '__main__':
    main()
