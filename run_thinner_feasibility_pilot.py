# -*- coding: utf-8 -*-
"""Research-feasibility pilot: real ChemExpo formulation + hypothetical benzene trace.
All numbers are exploratory ideal-liquid equilibrium calculations, NOT observations.
Requires pandas, numpy. No AI, no fitted data, no claimed worker exposure.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import math
D=Path(__file__).resolve().parent
# ChemExpo 04 입력 파일은 사용자가 지정한 Windows 폴더에서 읽습니다.
INPUT_DIR = Path(r'D:\fire\chemexpo_paint_thinner_output')
inp = INPUT_DIR / '04_model_input_candidates_REVIEW.csv'
if not inp.is_file():
    raise FileNotFoundError(f'ChemExpo 04 입력 파일을 찾을 수 없습니다: {inp}')
df=pd.read_csv(inp, dtype={'CAS':str,'Data Document ID':str})
# Confirmed NIST Antoine coefficients bar, Kelvin with valid T-range; EPA Perry coefficients Pa, Kelvin
props={
 '71-43-2':('벤젠','Benzene',78.1118,0.1,'NIST_Antoine',(4.01814,1203.835,-53.226),
            'https://webbook.nist.gov/cgi/cbook.cgi?ID=C71432&Type=ANTOINE', 'https://www.cdc.gov/niosh/npg/npgd0049.html'),
 '67-64-1':('아세톤','Acetone',58.0791,250.,'NIST_Antoine',(4.42448,1312.253,-32.445),
            'https://webbook.nist.gov/cgi/cbook.cgi?ID=C67641&Type=ANTOINE','https://www.cdc.gov/niosh/npg/npgd0004.html'),
 '141-78-6':('에틸아세테이트','Ethyl acetate',88.1051,400.,'NIST_Antoine',(4.22809,1245.702,-55.189),
            'https://webbook.nist.gov/cgi/cbook.cgi?ID=C141786&Type=ANTOINE','https://www.cdc.gov/niosh/npg/npgd0260.html'),
 '108-88-3':('톨루엔','Toluene',92.1384,100.,'NIST_Antoine',(4.14157,1377.578,-50.507),
            'https://webbook.nist.gov/cgi/cbook.cgi?ID=C108883&Type=ANTOINE','https://www.cdc.gov/niosh/npg/npgd0619.html'),
 '67-56-1':('메탄올','Methanol',32.0419,200.,'NIST_Antoine',(5.20409,1581.341,-33.50),
            'https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Type=ANTOINE','https://www.cdc.gov/niosh/npg/npgd0397.html'),
 '123-86-4':('n-부틸아세테이트','n-Butyl acetate',116.158,150.,'EPA_Perry',(71.34,-7285.8,-6.9459,9.99e-18,6),
            'https://19january2021snapshot.epa.gov/sites/static/files/documents/Romic_Section_D.pdf','https://www.cdc.gov/niosh/npg/npgd0072.html'),
 '108-10-1':('메틸이소부틸케톤','MIBK',100.1589,50.,'EPA_Perry',(153.23,-10055.,-19.848,1.64e-05,2),
            'https://19january2021snapshot.epa.gov/sites/static/files/documents/Romic_Section_D.pdf','https://www.cdc.gov/niosh/npg/npgd0326.html'),
}
# 3 formulations where EPA entries are distinct chemicals, sum 100; use actual extracted rows
ids=['1208563','273605','1183746']
labels={'1208563':'EPA-1 / Silver paint thinner','273605':'EPA-2 / Lacquer thinner','1183746':'EPA-3 / Vinyl thinner'}
T=298.15
C_total=100. # hypothetical ppm, same for all, not measured
benzene_observed=[0.039,0.040] # from separate 2022 paper; not measured in EPA products
rows=[]; detailed=[]
def vpkpa(cas):
 k,name,mw,oel,mode,par,*_=props[cas]
 if mode=='NIST_Antoine':
  a,b,c=par
  return 100*10**(a-b/(T+c)) # bar to kPa
 a,b,c,d,e=par
 return math.exp(a+b/T+c*math.log(T)+d*T**e)/1000

for pid in ids:
 base=df[df['Data Document ID']==pid].copy()
 assert len(base)>=2 and set(base.CAS).issubset(props)
 assert abs(base.wt_percent.sum()-100)<1e-8
 for dose in [0.,0.039,0.040,0.099,0.100]:
  comps=[(r.CAS, float(r.wt_percent)*(100-dose)/100.0,'EPA_SDS_reported') for _,r in base.iterrows()]
  if dose: comps.append(('71-43-2',dose,'2022_measured_value_transferred_as_counterfactual' if dose in benzene_observed else 'hypothetical_threshold_sensitivity'))
  ns=[wt/props[cas][2] for cas,wt,_ in comps]
  xs=np.array(ns)/sum(ns)
  partial=np.array([x*vpkpa(cas) for x,(cas,wt,source) in zip(xs,comps)])
  ys=partial/partial.sum()
  rates=np.array([y*1e6/props[cas][3] for y,(cas,wt,source) in zip(ys,comps)]) # vapor normalized per 1e6 ppm; relative only
  contribution=rates/rates.sum()
  for i,(cas,wt,source) in enumerate(comps):
   k,name,mw,oel,mode,par,vsrc,osrc=props[cas]
   detail={'product_id':pid,'product_label':labels[pid], 'scenario_benzene_wt_pct':dose,
        'cas':cas,'chemical':name,'liquid_wt_pct':wt,'liquid_mol_fraction':xs[i],
        'pure_vapor_pressure_25C_kPa':vpkpa(cas),'partial_vapor_pressure_kPa':partial[i],
        'vapor_mol_pct':ys[i]*100,'vapor_enrichment_vs_liquid_wt_ratio':ys[i]*100/wt,
        'oel_NIOSH_REL_ppm':oel,'hypothetical_air_ppm_at_100ppm_total':C_total*ys[i],
        'normalized_hazard_weight':rates[i], 'relative_contribution_pct':contribution[i]*100,
        'liquid_below_0p1_pct': bool(wt<0.1),
        'component_basis':source,'VP_source':vsrc,'OEL_source':osrc}
   detailed.append(detail)
  if dose>0:
   bz_i=[j for j,com in enumerate(comps) if com[0]=='71-43-2'][0]
   alt_rates=rates.copy();alt_rates[bz_i]/=10  # OSHA benzene PEL 1 ppm instead of NIOSH REL 0.1 ppm; other OELs unchanged
   rows.append({'product_id':pid,'product_label':labels[pid], 'EPA_base_components':len(base),
     'injected_benzene_wt_pct':dose, 'benzene_under_0p1pct':dose<0.1,
     'benzene_liquid_mol_pct':xs[bz_i]*100,'benzene_vapor_mol_pct':ys[bz_i]*100,
     'benzene_vapor_over_liquid_wt_fold':ys[bz_i]*100/dose,
     'benzene_relative_hazard_contribution_NIOSH_pct':contribution[bz_i]*100,
     'benzene_relative_hazard_contribution_OSHA_benzene1ppm_sensitivity_pct':alt_rates[bz_i]/alt_rates.sum()*100,
     'assumed_total_air_vapor_ppm':C_total,'hypothetical_benzene_air_ppm':ys[bz_i]*C_total,
     'hypothetical_Bz_over_NIOSH_REL':ys[bz_i]*C_total/0.1,
     'basis': 'EPA product composition + independent 2022 benzene concentration in counterfactual, not measured mixture',
     'assumptions':'Raoult ideal liquid gamma=1, equilibrium headspace 25C, constant 100ppm hypothetical airborne total; no field validation'
    })
summary=pd.DataFrame(rows)
detail=pd.DataFrame(detailed)
summary.to_csv(D/'thinner_pilot_summary.csv',index=False,encoding='utf-8-sig')
detail.to_csv(D/'thinner_pilot_component_detail.csv',index=False,encoding='utf-8-sig')
refrows=[]
for cas,(kor,en,mw,oel,mode,pars,vsrc,osrc) in props.items():
 refrows.append({'CAS':cas,'chemical':kor,'Molecular_weight_g_mol':mw,'VP_25C_kPa':vpkpa(cas), 'NIOSH_REL_TWA_ppm':oel,
 'VP_source_type':mode,'VP_source_URL':vsrc,'NIOSH_REL_source_URL':osrc,'parameter_tuple':str(pars)})
pd.DataFrame(refrows).to_csv(D/'thinner_pilot_property_sources.csv',index=False,encoding='utf-8-sig')
print(summary[['product_id','injected_benzene_wt_pct','benzene_vapor_mol_pct','benzene_relative_hazard_contribution_NIOSH_pct','benzene_relative_hazard_contribution_OSHA_benzene1ppm_sensitivity_pct','hypothetical_benzene_air_ppm']].to_string(index=False,float_format=lambda x:f'{x:.4f}'))
print('VAPOR BASE (0 dose)')
print(detail[detail.scenario_benzene_wt_pct==0][['product_id','chemical','liquid_wt_pct','vapor_mol_pct','relative_contribution_pct']].to_string(index=False,float_format=lambda x:f'{x:.3f}'))
# Counterfactual activity-coefficient sensitivity for benzene only (not a validated gamma model)
sens=[]
for pid in ids:
 for dose in [0.039,0.040,0.099]:
  base=df[df['Data Document ID']==pid]
  cs=[(r.CAS,float(r.wt_percent)*(100-dose)/100) for _,r in base.iterrows()]+[('71-43-2',dose)]
  ns=np.array([wt/props[cas][2] for cas,wt in cs]); xs=ns/ns.sum()
  for gamma in [0.5,1.,2.]:
   pp=np.array([xs[i]*vpkpa(cas)*(gamma if cas=='71-43-2' else 1.) for i,(cas,_) in enumerate(cs)])
   ys=pp/pp.sum()
   for oel in [0.1,1.0]:
    hazard=np.array([ys[i]/(oel if cas=='71-43-2' else props[cas][3]) for i,(cas,_) in enumerate(cs)])
    bzidx=len(cs)-1
    sens.append({'product_id':pid,'product_label':labels[pid],'benzene_wt_pct':dose,'gamma_benzene_scenario':gamma,
    'benzene_reference_limit_ppm':oel,'benzene_limit_basis':'NIOSH_REL' if oel==0.1 else 'OSHA_PEL_benzene_only',
    'benzene_vapor_mol_pct':ys[bzidx]*100,'benzene_contribution_pct':hazard[bzidx]/hazard.sum()*100,
    'scenario_not_measured':True})
pd.DataFrame(sens).to_csv(D/'thinner_pilot_gamma_sensitivity.csv',index=False,encoding='utf-8-sig')
print('BENZENE GAMMA SENSITIVITY @0.039 & NIOSH REL:')
print(pd.DataFrame(sens).query('benzene_wt_pct==0.039 and benzene_reference_limit_ppm==0.1')[['product_id','gamma_benzene_scenario','benzene_vapor_mol_pct','benzene_contribution_pct']].to_string(index=False,float_format=lambda x:f'{x:.3f}'))
