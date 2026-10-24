# -*- coding: utf-8 -*-
"""Finite liquid supply robustness test for illustrative well-mixed-room scenarios.

The existing analytical source assumes replenished constant-composition liquid
throughout the stated activity time. This audit computes the *implied*
evaporation mass using modeled species equilibrium vapor density and tests
whether a hypothetical 1 kg or 5 kg supply would suffice.

When exhausted, the source is switched off at the implied exhaustion time
and the room air decays under ventilation. The composition-invariance,
constant mass-transfer coefficient, no aerosol, and perfect-mixing assumptions
remain unvalidated and MUST NOT be described as measured workplace exposure.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import analyze_v2 as a2
from analyze_ventilated_scenarios import (
    SETTINGS, ROOM_M3, OBS_H, factors, OUT
)

R_GAS = 8.314462618  # Pa m3/(mol K)
MASSES_KG = (1.0, 5.0)


def mass_required_kg(t_h, area, k_m_h, ach_h, vapour_kg_per_m3,
                     volume_m3=ROOM_M3):
    """Integrated kg evaporated while replenished film is present."""
    if min(t_h, area, k_m_h, ach_h, vapour_kg_per_m3) < 0:
        raise ValueError('nonnegative parameters required')
    beta = area * k_m_h / volume_m3
    rate = beta + ach_h
    if t_h <= 0 or beta <= 0 or rate <= 0:
        return 0.0
    css = beta / rate
    integral_air_fraction = css * (
        t_h - (-np.expm1(-rate*t_h)) / rate)
    return float(area * k_m_h * vapour_kg_per_m3 *
                 (t_h - integral_air_fraction))


def truncate_to_stock_kg(duration_h, stock_kg, area, k_m_h, ach_h,
                         vapour_kg_per_m3):
    """Find source termination due to stock depletion (monotonic inverse)."""
    if stock_kg < 0:
        raise ValueError('Stock must be nonnegative')
    available = mass_required_kg(
        duration_h, area, k_m_h, ach_h, vapour_kg_per_m3)
    if stock_kg >= available:
        return float(duration_h), available, False
    low, high = 0., float(duration_h)
    for _ in range(60):
        mid = (low + high)/2
        used = mass_required_kg(mid, area, k_m_h, ach_h,
                                vapour_kg_per_m3)
        if used > stock_kg:
            high = mid
        else:
            low = mid
    return (low+high)/2, available, True


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(a2.RES/'species_central_270.csv.gz',
                     dtype={'CAS':str})
    sp = sp[sp['model'].eq('HANNA') & sp['T'].isin(('25C','35C'))].copy()
    props = pd.read_csv(a2.V2/'data'/'pure_properties.csv',
                        encoding='utf-8-sig',dtype={'CAS':str})
    mw = props.set_index('CAS')['MW_g_mol']
    sp['MW_g_mol'] = sp['CAS'].map(mw)
    if sp['MW_g_mol'].isna().any():
        raise ValueError('No mass density for unresolved modeled species')
    sp['temperature_K'] = sp['T'].map({'25C':298.15,'35C':308.15})
    sp['vapor_density_kg_per_m3'] = (
        sp['p_kPa'] * 1000 * sp['MW_g_mol'] /
        (R_GAS * sp['temperature_K'] * 1000))
    density = sp.groupby(['문서키','T'],as_index=False)[
        'vapor_density_kg_per_m3'].sum()
    baseline = pd.read_csv(OUT/'V1_wellmixed_270_by_scenario.csv',
                           encoding='utf-8-sig')
    if baseline['문서키'].nunique()!=270:
        raise ValueError('Missing room model on 270 documents')
    d=baseline.merge(density,left_on=['문서키','temperature'],
                     right_on=['문서키','T'],validate='many_to_one')
    if len(d)!=len(baseline):
        raise ValueError('Missing molecular-vapor density')
    params = {x[0]:x[1:] for x in SETTINGS}
    rows=[]
    for row in d.itertuples(index=False):
        area,k,ach,requested_h=params[row.scenario_name]
        density_k=float(row.vapor_density_kg_per_m3)
        for stock in MASSES_KG:
            actual_h, required_kg, exhausted=truncate_to_stock_kg(
                requested_h,stock,area,k,ach,density_k)
            _,f8=factors(area,k,ach,actual_h)
            rows.append({
                '문서키':row.문서키,
                '기존제품ID':row.기존제품ID,
                '조성지문':row.조성지문,
                'temperature':row.temperature,
                'scenario_name':row.scenario_name,
                'hypothetical_inventory_kg':stock,
                'required_mass_if_replenished_kg':required_kg,
                'inventory_exhausted_before_source_period_end':exhausted,
                'effective_source_duration_h':actual_h,
                'intended_source_duration_h':requested_h,
                'equilibrium_vapor_density_kg_per_m3':density_k,
                'finite_mass_mean_fraction_of_saturation':f8,
                'finite_mass_hypo8h_TierA_index':row.I_sat_TierA*f8,
                'finite_mass_hypo8h_TierA_index_P5':row.I_sat_TierA_P5*f8,
                'OELcov_TierA':row.OELcov_TierA,
                'actual_workplace_exposure_measured':False,
                'formula_limitation':'assumes unchanged solvent composition until exhausted'
            })
    result=pd.DataFrame(rows)
    result.to_csv(OUT/'V5_finite_source_mass_sensitivity.csv',
                  index=False,encoding='utf-8-sig')
    summary=[]
    for (temp,sc,stock),x in result.groupby([
        'temperature','scenario_name','hypothetical_inventory_kg']):
        eligible=x[x['OELcov_TierA']>=0.9]
        summary.append({
            'temperature':temp,'scenario_name':sc,
            'hypothetical_inventory_kg':stock,
            'documents_n':len(x),'OEL_90pct_eligible_n':len(eligible),
            'source_depleted_documents_n':int(x[
                'inventory_exhausted_before_source_period_end'].sum()),
            'finite_mass_hypo8h_index_ge1_n':int(eligible[
                'finite_mass_hypo8h_TierA_index'].ge(1).sum()),
            'finite_mass_P5_ge1_n':int(eligible[
                'finite_mass_hypo8h_TierA_index_P5'].ge(1).sum()),
            'finite_mass_hypo8h_index_ge10_n':int(eligible[
                'finite_mass_hypo8h_TierA_index'].ge(10).sum()),
            'median_required_mass_if_replenished_kg':float(
                x['required_mass_if_replenished_kg'].median()),
            'is_experimental_validation':False
        })
    pd.DataFrame(summary).to_csv(
        OUT/'V6_mass_budget_scenario_summary.csv',index=False,
        encoding='utf-8-sig')
    headline={
        'basis':'analytic room mass balance plus finite-amount supply sensitivity',
        'stocks_kg':list(MASSES_KG),
        'source_temperature_assumption':'constant 25C/35C vapor density during source period',
        'reminder':'hypothetical not observed, no finite-composition evolution',
        'summary':summary
    }
    (OUT/'headline_mass_budget.json').write_text(
        json.dumps(headline,ensure_ascii=False,indent=2),encoding='utf-8')
    print('HYPOTHETICAL_FINITE_SOURCE_MASS '+
          json.dumps(headline,ensure_ascii=False))


if __name__=='__main__':
    main()
