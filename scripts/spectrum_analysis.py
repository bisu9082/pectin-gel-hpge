# -*- coding: utf-8 -*-
"""pectin JER revision — spectrum reanalysis of 26 CNF files.

For each spectrum and each gamma line:
  net peak area (side-band linear background), counting uncertainty,
  net count rate (per live time), Currie critical level & detection limit.
Then paired recovery per line (efficiency-free rate ratio, with mass and
decay corrections) and multi-line consistency for Cs-134 / Co-60.
"""
import sys, os, glob, math, json, datetime
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnf_reader import read_cnf

BASE = '/tmp/pectin_new/펙틴 스펙트럼'
OUT = '/home/claude/pectin_rev/out'
os.makedirs(OUT, exist_ok=True)

LINES = [  # nuclide, energy keV, branching
    ('Am-241',   59.54, 0.359),
    ('Cs-134',  604.72, 0.976),
    ('Cs-134',  795.86, 0.855),
    ('Cs-137',  661.66, 0.851),
    ('Co-60',  1173.23, 0.9985),
    ('Co-60',  1332.49, 0.9998),
]
HALF_LIFE_Y = {'Cs-137': 30.08, 'Cs-134': 2.0652, 'Am-241': 432.6, 'Co-60': 5.2711}

# sample masses (g) from 펙틴 실험_0824.xlsx
MASS = {
    ('orig', 1): (28.2664, 28.3941), ('orig', 2): (28.2692, 28.7904),
    ('orig', 3): (28.2665, 29.1465), ('orig', 4): (28.2675, 29.2001),
    ('orig', 5): (28.2664, 29.0141),
    ('repro', 1): (28.2685, 29.2117), ('repro', 2): (28.2683, 29.0428),
    ('repro', 3): (28.2665, 29.2771),
}

FILES = {}
for i in range(1, 6):
    FILES[('orig', i, 'before')] = os.path.join(BASE, '고형화전', f'{i}번 시료_25년 보정.CNF')
    FILES[('orig', i, 'after')]  = os.path.join(BASE, '고형화후', f'{i}번 시료_분석완료.CNF')
for i in range(1, 4):
    FILES[('repro', i, 'before')] = os.path.join(BASE, '재현성테스트', '고형화 전', f'{i}.CNF')
    FILES[('repro', i, 'after')]  = os.path.join(BASE, '재현성테스트', '고형화 후', f'{i}.CNF')


def ch_of_energy(E, cc, nch):
    """Invert (possibly cubic) energy calibration numerically."""
    ch = np.arange(nch)
    en = cc[0] + cc[1]*ch + cc[2]*ch**2 + cc[3]*ch**3
    return float(np.interp(E, en, ch))


def peak_net(counts, cc, E, pk_hw=3.5, bg_gap=4.5, bg_w=4.5):
    """Net area via side-band background subtraction.
    Windows in keV: peak E±pk_hw; bg [E-bg_gap-bg_w, E-bg_gap] and mirror."""
    nch = len(counts)
    c = lambda e: int(round(ch_of_energy(e, cc, nch)))
    p1, p2 = c(E-pk_hw), c(E+pk_hw)
    l1, l2 = c(E-bg_gap-bg_w), c(E-bg_gap)
    r1, r2 = c(E+bg_gap), c(E+bg_gap+bg_w)
    G = counts[p1:p2+1].sum()
    n_pk = p2 - p1 + 1
    BL, BR = counts[l1:l2+1].sum(), counts[r1:r2+1].sum()
    n_bg = (l2-l1+1) + (r2-r1+1)
    scale = n_pk / n_bg
    B = scale * (BL + BR)
    net = G - B
    u_net = math.sqrt(G + B*scale)          # var(G)+var(B_est)
    # Currie: critical level and detection limit (counts)
    Lc = 2.326 * math.sqrt(2*B*scale) if B > 0 else 2.326
    Ld = 2.71 + 4.65*math.sqrt(B*scale) if B > 0 else 2.71
    return dict(gross=float(G), bg=float(B), net=float(net), u_net=float(u_net),
                Lc=float(Lc), Ld=float(Ld), n_pk=n_pk)


def main():
    rows = []
    spec = {}
    for key, fn in FILES.items():
        d = read_cnf(fn)
        spec[key] = d
        for nuc, E, br in LINES:
            r = peak_net(d['counts'], d['ecal'], E)
            rate = r['net'] / d['livetime']
            u_rate = r['u_net'] / d['livetime']
            rows.append(dict(group=key[0], sample=key[1], phase=key[2],
                             nuclide=nuc, E=E, br=br,
                             start=d['start'].isoformat(), lt=d['livetime'],
                             **r, rate=rate, u_rate=u_rate,
                             detected=r['net'] > r['Lc']))
    import csv
    with open(f'{OUT}/peaks_all.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # paired recovery per line (activity-concentration basis, Bq/kg-of-sample-mass)
    def get(group, s, phase, nuc, E):
        for r in rows:
            if (r['group'], r['sample'], r['phase'], r['nuclide'], r['E']) == (group, s, phase, nuc, E):
                return r
    rec = []
    for (group, s), (mb, ma) in MASS.items():
        tb = spec[(group, s, 'before')]['start']
        ta = spec[(group, s, 'after')]['start']
        dt_y = (ta - tb).total_seconds() / (365.25*86400)
        for nuc, E, br in LINES:
            b = get(group, s, 'before', nuc, E)
            a = get(group, s, 'after', nuc, E)
            if b['net'] <= b['Lc'] or a['net'] <= a['Lc']:
                rec.append(dict(group=group, sample=s, nuclide=nuc, E=E,
                                R=None, uR=None, note='below Lc (b:%d a:%d)' % (b['detected'], a['detected'])))
                continue
            decay = math.exp(math.log(2)/HALF_LIFE_Y[nuc]*dt_y)  # correct after back to before date
            # two conventions:
            R_rate = a['rate']/b['rate']*decay*100.0                      # same-mass basis
            R_conc = a['rate']/ma/(b['rate']/mb)*decay*100.0              # per-mass basis (Genie Bq/kg)
            uR = R_conc*math.sqrt((a['u_rate']/a['rate'])**2+(b['u_rate']/b['rate'])**2)
            rec.append(dict(group=group, sample=s, nuclide=nuc, E=E,
                            R_rate=round(R_rate,2), R_conc=round(R_conc,2), uR=round(uR,2),
                            dt_d=round(dt_y*365.25,1), note=''))
    with open(f'{OUT}/recovery_lines.csv', 'w', newline='') as f:
        keys = ['group','sample','nuclide','E','R_rate','R_conc','uR','dt_d','note','R','uR2']
        w = csv.DictWriter(f, fieldnames=['group','sample','nuclide','E','R_rate','R_conc','uR','dt_d','note'], extrasaction='ignore')
        w.writeheader(); w.writerows(rec)

    # print summaries
    print('=== 회수율 (라인별, R_conc = per-mass basis / R_rate = same-mass) ===')
    for r in rec:
        if r.get('R_conc') is not None:
            print(f"{r['group']:5s} S{r['sample']} {r['nuclide']:7s} {r['E']:7.1f} keV  "
                  f"R_conc={r['R_conc']:7.2f}±{r['uR']:5.2f}%  R_rate={r['R_rate']:7.2f}%")
        else:
            print(f"{r['group']:5s} S{r['sample']} {r['nuclide']:7s} {r['E']:7.1f} keV  -- {r['note']}")
    print()
    print('=== 검출 요약 (net, Lc, 검출여부) 후단 확인용 ===')
    for r in rows:
        if not r['detected']:
            print(f"UNDET {r['group']:5s} S{r['sample']} {r['phase']:6s} {r['nuclide']:7s} {r['E']:7.1f}  net={r['net']:7.1f} Lc={r['Lc']:6.1f}")

if __name__ == '__main__':
    main()
