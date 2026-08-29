# -*- coding: utf-8 -*-
"""Definitive dataset rebuild from 펙틴 실험_0824.xlsx raw Genie values.
All values transcribed from the xlsx (시료별정리 + 재현성 sheets).
Computes: paired recovery R ± u (k=1), |Δ| vs k=2 test, per-nuclide stats,
activity-level comparison, MDA summary. Outputs CSV + console summary."""
import numpy as np, math, csv
from scipy import stats as st

# (nuclide, set, sample, A_KINS_Cs137basis?, A_before, u_b, MDA_b, A_after, u_a, MDA_a)
# KINS certified reference values apply to the seawater matrix per round:
# round 2024: 45.1±0.7 (S1,S2), round 2022: 24.6±0.3 (S3,S4, repro1-3), round 2019: 24.8±0.2 (S5)
# NOTE: those certified values are for Cs-137 [확인: 타핵종 기준값은 KINS 성적서 필요]
D = [
 # nuclide, grp, s, Ab, ub, MDAb, Aa, ua, MDAa
 ('Cs-137','orig',1, 43.885,1.507,2.916, 44.279,1.533,3.246),
 ('Cs-137','orig',2, 47.119,1.552,2.444, 43.286,1.526,3.053),
 ('Cs-137','orig',3, 30.527,1.314,2.464, 30.474,1.299,2.650),
 ('Cs-137','orig',4, 28.035,1.134,4.312, 30.280,1.277,2.999),
 ('Cs-137','orig',5, 30.101,1.295,2.570, 29.684,2.980,3.140),
 ('Cs-137','repro',1, 26.360,1.222,2.92, 24.210,1.177,3.00),
 ('Cs-137','repro',2, 28.420,1.260,2.51, 27.710,1.267,None),
 ('Cs-137','repro',3, 27.770,0.784,1.705, 28.150,1.258,2.720),
 ('Cs-134','orig',1, 43.287,1.319,3.553, 39.554,1.337,3.985),
 ('Cs-134','orig',2, 41.811,1.326,3.570, 42.053,1.369,4.653),
 ('Cs-134','orig',3, 27.799,1.643,7.270, 25.645,1.557,5.110),
 ('Cs-134','orig',4, 29.610,1.645,2.285, 24.170,1.542,5.362),
 ('Cs-134','repro',1, 27.340,1.785,7.190, 25.430,1.805,7.984),
 ('Cs-134','repro',2, 25.820,1.736,7.402, 26.510,1.810,6.791),
 ('Cs-134','repro',3, 26.620,1.121,4.272, 26.380,1.830,7.126),
 ('Am-241','orig',1, 97.301,2.917,5.124, 93.185,2.895,5.333),
 ('Am-241','orig',2, 97.135,2.881,5.091, 96.385,2.909,5.226),
 ('Am-241','orig',3, 58.504,2.113,4.556, 50.354,2.049,4.280),
 ('Am-241','orig',4, 57.147,1.892,4.084, 52.681,2.070,4.606),
 ('Am-241','orig',5, 11.613,1.118,4.025, 11.683,1.189,7.637),
 ('Am-241','repro',1, 49.120,1.940,4.138, 50.200,2.049,4.555),
 ('Am-241','repro',2, 48.650,1.990,4.434, 46.200,1.948,4.210),
 ('Am-241','repro',3, 51.400,1.430,2.599, 48.580,1.970,4.259),
 ('Co-60','orig',1, 53.908,1.528,3.200, 51.024,1.503,3.796),
 ('Co-60','orig',2, 56.847,1.543,3.370, 53.970,1.528,2.801),
 ('Co-60','orig',5, 19.514,1.350,5.447, 18.990,1.828,5.883),
]
LEVEL = {('orig',1):'high45',('orig',2):'high45',('orig',3):'low25',('orig',4):'low25',
         ('orig',5):'low25',('repro',1):'low25',('repro',2):'low25',('repro',3):'low25'}

rows=[]
for nuc,grp,s,Ab,ub,MDAb,Aa,ua,MDAa in D:
    R = Aa/Ab*100
    uR = R*math.sqrt((ua/Aa)**2+(ub/Ab)**2)
    delta = R-100
    within_k2 = abs(delta) <= 2*uR
    rows.append(dict(nuclide=nuc,set=grp,sample=s,level=LEVEL[(grp,s)],
                     A_before=Ab,u_before=ub,MDA_before=MDAb,
                     A_after=Aa,u_after=ua,MDA_after=MDAa,
                     R=round(R,1),uR=round(uR,1),Delta=round(delta,1),
                     consistent_k2=within_k2, in_80_120=(80<=R<=120)))

with open('/home/claude/pectin_rev/out/dataset_definitive.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print('=== 26쌍 확정 데이터셋 (Genie 원자료 기준) ===')
for r in rows:
    tag = 'OK ' if r['consistent_k2'] else '**'
    print(f"{r['nuclide']:7s} {r['set']:5s} S{r['sample']} [{r['level']}] "
          f"{r['A_before']:7.3f}±{r['u_before']:.3f} -> {r['A_after']:7.3f}±{r['u_after']:.3f}  "
          f"R={r['R']:6.1f}±{r['uR']:4.1f}%  Δ={r['Delta']:+6.1f}  k2:{tag} 80-120:{r['in_80_120']}")

print()
print('=== 핵종별 통계 ===')
summ={}
for nuc in ['Cs-137','Cs-134','Am-241','Co-60']:
    v=np.array([r['R'] for r in rows if r['nuclide']==nuc])
    n=len(v); m=v.mean(); sd=v.std(ddof=1)
    tcrit=st.t.ppf(0.975,n-1); ci=(m-tcrit*sd/math.sqrt(n), m+tcrit*sd/math.sqrt(n))
    summ[nuc]=(n,m,sd,ci,v.min(),v.max())
    print(f'{nuc:7s}: n={n}  mean={m:6.2f} ± {sd:5.2f} (SD)  95%CI=[{ci[0]:.1f},{ci[1]:.1f}]  min={v.min():.1f} max={v.max():.1f}')
allR=np.array([r['R'] for r in rows])
print(f'전체 26쌍: 80–120% 내 {sum(1 for r in rows if r["in_80_120"])}/26, |Δ|<=2u(k=2 일치) {sum(1 for r in rows if r["consistent_k2"])}/26')

print()
print('=== 활동도 레벨 비교 (~45 vs ~25 Bq/kg, 전핵종 풀링) ===')
hi=np.array([r['R'] for r in rows if r['level']=='high45'])
lo=np.array([r['R'] for r in rows if r['level']=='low25'])
t,p = st.ttest_ind(hi,lo,equal_var=False)
u,pu = st.mannwhitneyu(hi,lo,alternative='two-sided')
print(f'high45: n={len(hi)} mean={hi.mean():.2f}±{hi.std(ddof=1):.2f} | low25: n={len(lo)} mean={lo.mean():.2f}±{lo.std(ddof=1):.2f}')
print(f'Welch t={t:.3f} p={p:.3f} | Mann-Whitney U={u:.0f} p={pu:.3f} | 차이={hi.mean()-lo.mean():+.2f}%p')

print()
print('=== MDA 요약 (Genie, Bq/kg) ===')
for nuc in ['Cs-137','Cs-134','Am-241','Co-60']:
    mb=[r['MDA_before'] for r in rows if r['nuclide']==nuc and r['MDA_before']]
    ma=[r['MDA_after'] for r in rows if r['nuclide']==nuc and r['MDA_after']]
    print(f'{nuc:7s}: before {min(mb):.1f}–{max(mb):.1f}  after {min(ma):.1f}–{max(ma):.1f}')

print()
print('=== 원고 대비 변경점 ===')
old={'Cs-137':(5,99.9,5.6),'Cs-134':(4,88.8,12.6),'Am-241':(5,93.7,5.6),'Co-60':(3,94.8,2.7)}
for nuc,(n0,m0,s0) in old.items():
    n,m,sd,ci,mn,mx = summ[nuc]
    print(f'{nuc:7s}: 원고 n={n0} {m0}±{s0} -> 신규 n={n} {m:.1f}±{sd:.1f}')
