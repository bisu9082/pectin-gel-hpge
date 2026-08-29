# -*- coding: utf-8 -*-
"""Definitive 26-pair dataset and the statistics reported in the manuscript.

Activity concentrations and their combined standard uncertainties (k=1) and
MDA values are transcribed from the Genie-2000 analysis outputs; the sample
masses are those entered in the same analyses. Paired retention is

    R (%) = (A_after * m_after) / (A_before * m_before) * 100

Outputs (raw_data/): dataset_retention.csv
Reproduces: Table 1, Table 2, Table 3 and the statistics of Sections 3.3-3.5.
"""
import csv, math
import numpy as np
from scipy import stats as st
from _paths import RAW

# nuclide, sample, A_before, u_before, MDA_before, A_after, u_after, MDA_after
D = [
 ('Cs-137', 1, 43.885, 1.507, 2.916, 44.279, 1.533, 3.246),
 ('Cs-137', 2, 47.119, 1.552, 2.444, 43.286, 1.526, 3.053),
 ('Cs-137', 3, 30.527, 1.314, 2.464, 30.474, 1.299, 2.650),
 ('Cs-137', 4, 28.035, 1.134, 4.312, 30.280, 1.277, 2.999),
 ('Cs-137', 5, 30.101, 1.295, 2.570, 29.684, 2.980, 3.140),
 ('Cs-137', 6, 26.360, 1.222, 2.920, 24.210, 1.177, 3.000),
 ('Cs-137', 7, 28.420, 1.260, 2.510, 27.710, 1.267, None),
 ('Cs-137', 8, 27.770, 0.784, 1.705, 28.150, 1.258, 2.720),
 ('Cs-134', 1, 43.287, 1.319, 3.553, 39.554, 1.337, 3.985),
 ('Cs-134', 2, 41.811, 1.326, 3.570, 42.053, 1.369, 4.653),
 ('Cs-134', 3, 27.799, 1.643, 7.270, 25.645, 1.557, 5.110),
 ('Cs-134', 4, 29.610, 1.645, 2.285, 24.170, 1.542, 5.362),
 ('Cs-134', 6, 27.340, 1.785, 7.190, 25.430, 1.805, 7.984),
 ('Cs-134', 7, 25.820, 1.736, 7.402, 26.510, 1.810, 6.791),
 ('Cs-134', 8, 26.620, 1.121, 4.272, 26.380, 1.830, 7.126),
 ('Am-241', 1, 97.301, 2.917, 5.124, 93.185, 2.895, 5.333),
 ('Am-241', 2, 97.135, 2.881, 5.091, 96.385, 2.909, 5.226),
 ('Am-241', 3, 58.504, 2.113, 4.556, 50.354, 2.049, 4.280),
 ('Am-241', 4, 57.147, 1.892, 4.084, 52.681, 2.070, 4.606),
 ('Am-241', 5, 11.613, 1.118, 4.025, 11.683, 1.189, 7.637),
 ('Am-241', 6, 49.120, 1.940, 4.138, 50.200, 2.049, 4.555),
 ('Am-241', 7, 48.650, 1.990, 4.434, 46.200, 1.948, 4.210),
 ('Am-241', 8, 51.400, 1.430, 2.599, 48.580, 1.970, 4.259),
 ('Co-60',  1, 53.908, 1.528, 3.200, 51.024, 1.503, 3.796),
 ('Co-60',  2, 56.847, 1.543, 3.370, 53.970, 1.528, 2.801),
 ('Co-60',  5, 19.514, 1.350, 5.447, 18.990, 1.828, 5.883),
]
MASS = {1: (28.2664, 28.3941), 2: (28.2692, 28.7904), 3: (28.2665, 29.1465),
        4: (28.2675, 29.2001), 5: (28.2664, 29.0141), 6: (28.2685, 28.2685),
        7: (28.2683, 28.2683), 8: (28.2665, 28.2665)}
# S1-S2: 2024 round (45.1 Bq/kg); S3,S4,S6-S8: 2022 round (24.6); S5: 2019 (24.8)
LEVEL = {1: 'high45', 2: 'high45', 3: 'low25', 4: 'low25', 5: 'low25',
         6: 'low25', 7: 'low25', 8: 'low25'}
SET = {s: ('original' if s <= 5 else 'replicate') for s in range(1, 9)}
EGAMMA = {'Cs-137': 661.7, 'Cs-134': 604.7, 'Am-241': 59.5, 'Co-60': 1173.2}


def main():
    rows = []
    for nuc, s, Ab, ub, mdab, Aa, ua, mdaa in D:
        mb, ma = MASS[s]
        R = (Aa*ma)/(Ab*mb)*100
        uR = R*math.sqrt((ua/Aa)**2 + (ub/Ab)**2)
        rows.append(dict(nuclide=nuc, sample=f'S{s}', set=SET[s], level=LEVEL[s],
                         E_keV=EGAMMA[nuc], A_before=Ab, u_before=ub,
                         MDA_before=mdab, A_after=Aa, u_after=ua, MDA_after=mdaa,
                         m_before=mb, m_after=ma, R_percent=round(R, 1),
                         u_R=round(uR, 1), Delta_percent=round(R-100, 1),
                         consistent_k2=abs(R-100) <= 2*uR,
                         within_80_120=80 <= R <= 120))
    with open(RAW/'dataset_retention.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    print('=== Table 2: paired retention by radionuclide ===')
    for nuc in ['Cs-137', 'Cs-134', 'Am-241', 'Co-60']:
        v = np.array([r['R_percent'] for r in rows if r['nuclide'] == nuc])
        n = len(v); sd = v.std(ddof=1)
        ci = st.t.interval(0.95, n-1, loc=v.mean(), scale=sd/math.sqrt(n))
        print(f'  {nuc:7s} n={n}  {v.mean():6.1f} +- {sd:4.1f}  '
              f'95% CI [{ci[0]:.1f}, {ci[1]:.1f}]  min {v.min():.1f}  max {v.max():.1f}')

    print('\n=== Table 3: MDA ranges (Bq/kg) ===')
    for nuc in ['Cs-137', 'Cs-134', 'Am-241', 'Co-60']:
        mb = [r['MDA_before'] for r in rows if r['nuclide'] == nuc and r['MDA_before']]
        ma = [r['MDA_after'] for r in rows if r['nuclide'] == nuc and r['MDA_after']]
        print(f'  {nuc:7s} {EGAMMA[nuc]:7.1f} keV  before {min(mb):.1f}-{max(mb):.1f}'
              f'   after {min(ma):.1f}-{max(ma):.1f}')

    n_ok = sum(r['within_80_120'] for r in rows)
    n_k2 = sum(r['consistent_k2'] for r in rows)
    print(f'\n=== Section 3.3 ===\n  within 80-120%: {n_ok}/{len(rows)}'
          f'\n  consistent with 100% at k=2: {n_k2}/{len(rows)}')
    print('  statistically resolved changes: ' + ', '.join(
        f"{r['nuclide']} {r['sample']} ({r['R_percent']}+-{r['u_R']}%)"
        for r in rows if not r['consistent_k2']))

    hi = np.array([r['R_percent'] for r in rows if r['level'] == 'high45'])
    lo = np.array([r['R_percent'] for r in rows if r['level'] == 'low25'])
    t, p = st.ttest_ind(hi, lo, equal_var=False)
    u, pu = st.mannwhitneyu(hi, lo, alternative='two-sided')
    print(f'\n=== Activity-level comparison ===\n  ~45 Bq/kg: {hi.mean():.1f} +- '
          f'{hi.std(ddof=1):.1f} (n={len(hi)})   ~25 Bq/kg: {lo.mean():.1f} +- '
          f'{lo.std(ddof=1):.1f} (n={len(lo)})\n  difference {hi.mean()-lo.mean():+.1f} '
          f'percentage points; Welch p={p:.2f}, Mann-Whitney p={pu:.2f}')

    co = np.array([r['R_percent'] for r in rows if r['nuclide'] == 'Co-60'])
    rng = np.random.default_rng(20260829)
    for shift in np.arange(4.0, 12.01, 0.1):
        x = rng.normal(100+shift, co.std(ddof=1), size=(200000, 3))
        tt = (x.mean(1)-100)/(x.std(1, ddof=1)/math.sqrt(3))
        if (np.abs(tt) > st.t.ppf(0.975, 2)).mean() >= 0.80:
            print(f'\n=== SI S2: Co-60 sensitivity ===\n  80% power at a shift of '
                  f'~{shift:.1f} percentage points (n=3, alpha=0.05)')
            break
    print(f'\nWritten: {RAW/"dataset_retention.csv"}')


if __name__ == '__main__':
    main()
