# -*- coding: utf-8 -*-
"""Non-experimental hypothetical room-air sensitivity, 25 C / 35 C.

NOT a validated exposure model or a workplace assessment.

For a hypothetical well-mixed room, a replenished liquid surface, a single
mass-transfer coefficient shared by all components, and zero starting
background concentration:
    dC/dt = (k*A/V)*(Csat-C) - ACH*C
where Csat is the existing model-derived equilibrium saturation concentration
(ppm) of a compound, A [m2] is exposed liquid area, k [m/h] is the illustrative
mass-transfer coefficient, and ACH [1/h] is air exchange.

During a finite source period tau, and for a complete hypothetical 8h
observation including the ventilated decay afterward, the model has the
analytic solution used below. No solvent depletion, droplet aerosol, mixing
heterogeneity, multiple sources, adsorption, personal breathing-zone effects,
respirators, changing composition, or thermal gradients are represented.

The numerical choices are *exploratory parameters*, not measurements of small
workshops or standardized operating conditions. TWA-labelled indices are
model-derived 8h time means under those hypothetical premises, not evidence
of regulatory exceedance.
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

import analyze_v2 as a2

OUT = a2.RES / 'regulatory_gap'
OBS_H = 8.0
ROOM_M3 = 100.0

# A [m2], k [m/h], ACH [h-1], surface active time [h], V=100m3 fixed
SETTINGS = [
    ('small_surface_high_vent_15min', 0.1, 0.1, 8.0, 0.25),
    ('small_surface_moderate_vent_1h', 0.5, 0.2, 5.0, 1.0),
    ('medium_surface_moderate_vent_1h', 1.0, 0.5, 3.0, 1.0),
    ('large_surface_low_vent_4h', 3.0, 1.0, 1.0, 4.0),
    ('large_surface_very_low_vent_4h', 5.0, 1.0, 0.5, 4.0),
]


def factors(area_m2, k_m_h, ach_h, source_h, volume_m3=ROOM_M3,
            observe_h=OBS_H):
    """(source-end concentration / Csat, 8h mean / Csat)."""
    if (min(area_m2, k_m_h, ach_h, source_h, volume_m3, observe_h) < 0
            or volume_m3 <= 0 or observe_h <= 0 or source_h > observe_h):
        raise ValueError('Invalid room or activity assumptions')
    a = area_m2 * k_m_h / volume_m3
    removal = ach_h
    rate = a + removal
    tau = source_h
    if rate == 0 or tau == 0:
        return 0.0, 0.0
    css = a / rate
    end = css * (-np.expm1(-rate * tau))
    # Integral during release phase 0..tau.
    integral_source = css * (tau - (-np.expm1(-rate*tau)) / rate)
    # Integral after source closes: no k*A term and decay by air exchange.
    remaining = observe_h - tau
    integral_decay = (end * (-np.expm1(-removal * remaining)) / removal
                      if removal > 0 else end * remaining)
    avg = (integral_source + integral_decay) / observe_h
    if not (0 <= avg <= 1+1e-12 and 0 <= end <= 1+1e-12):
        raise AssertionError('Improper saturation-fraction bounds')
    return float(end), float(avg)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(a2.RES / 'scenarios_270.csv.gz')
    if 'scenario' not in s or not {'25C', '35C'}.issubset(set(s['T'])):
        raise ValueError('Expected both temperatures and the scenario identifier')
    summ = a2.summarise(s)
    base = summ[summ['model'].eq('HANNA') &
                summ['T'].isin(['25C','35C'])].copy()
    if base.groupby(['문서키','T']).size().max() != 1:
        raise AssertionError('Non-unique document-temperature summaries')
    _, prod, _ = a2.load()
    base = base.merge(prod[['문서키', '기존제품ID', '제조사', '제품명', '조성지문',
                             'UVCB대리물질사용']], on='문서키', validate='many_to_one')
    if base['문서키'].nunique() != 270:
        raise AssertionError('Expected the complete 270 document census')

    crosswalk = OUT / 'R1_constituent_rule_vapor_crosswalk.csv'
    coverage = OUT / 'R3_270_document_coverage.csv'
    if not crosswalk.is_file() or not coverage.is_file():
        raise FileNotFoundError(
            'Run analyze_regulatory_vapor_gap.py before ventilated sensitivity')
    r1 = pd.read_csv(crosswalk, encoding='utf-8-sig')
    r3 = pd.read_csv(coverage, encoding='utf-8-sig')
    below = r1[
        r1['기상모델'].eq('HANNA') &
        r1['해당규정_함량비교'].eq('below_rule_threshold') &
        r1['기상성분추적신뢰'].eq('central_model_value') &
        r1['모형_개별성분_TierA_OEL지수'].notna()
    ].copy()
    valid_rules = set(below['문서키'])
    row_list, rule_rows, scenario_inputs = [], [], []
    for key, area, transfer, ach, duration in SETTINGS:
        end, avg = factors(area, transfer, ach, duration)
        scenario_inputs.append({
            'scenario_name': key, 'room_volume_m3': ROOM_M3,
            'source_surface_m2': area, 'mass_transfer_m_per_h': transfer,
            'air_exchange_per_h': ach, 'source_duration_h': duration,
            'observation_h': OBS_H, 'end_of_use_to_saturation': end,
            'eight_hour_mean_to_saturation': avg,
            'assumption_scope': 'well-mixed, constant composition, replenished liquid, zero initial background'})
        for _, r in base.iterrows():
            row = {
                '문서키': r['문서키'], '기존제품ID': r['기존제품ID'],
                '제품명': r['제품명'], '조성지문': r['조성지문'],
                'temperature': r['T'], 'scenario_name': key,
                'OELcov_TierA': r['OELcov_TierA'],
                'OELcov_TierB': r['OELcov_TierB'],
                'I_sat_TierA': r['I_mix_TierA'],
                'I_sat_TierA_P5': r['I_mix_TierA_P5'],
                'I_sat_TierA_P95': r['I_mix_TierA_P95'],
                'I_sat_TierB': r['I_mix_TierB'],
                'end_use_saturation_fraction': end,
                'eight_hour_saturation_fraction': avg,
                'hypo_8h_TierA_index': r['I_mix_TierA']*avg,
                'hypo_8h_TierA_index_P5': r['I_mix_TierA_P5']*avg,
                'hypo_8h_TierA_index_P95': r['I_mix_TierA_P95']*avg,
                'hypo_8h_TierB_index': r['I_mix_TierB']*avg,
                'official_TWA_coverage_90pct': r['OELcov_TierA'] >= 0.9,
                'verified_rule_below_and_traceable': r['문서키'] in valid_rules,
                'source_note': 'No personal exposure, no validated real workplace air concentration'
            }
            row_list.append(row)
        # For selected CAS clauses, the 25C saturated single-component
        # contribution times hypothetical 8h dilution factor.
        tab = below[[
            '문서키', 'CAS_원문', '규정범주',
            'MSDS_액상하한_mass_pct', 'MSDS_액상상한_mass_pct',
            '규정_액상함량기준_mass_pct',
            '모형_개별성분_TierA_OEL지수']].copy()
        tab['scenario_name'] = key
        tab['temperature'] = '25C'
        tab['eight_hour_saturation_fraction'] = avg
        tab['hypo_8h_single_CAS_TierA_index'] = (
            tab['모형_개별성분_TierA_OEL지수'] * avg)
        rule_rows.append(tab)
    frames = pd.DataFrame(row_list)
    frames.to_csv(OUT / 'V1_wellmixed_270_by_scenario.csv',
                  index=False, encoding='utf-8-sig')
    pd.DataFrame(scenario_inputs).to_csv(
        OUT / 'V0_exposure_scenario_assumptions.csv',
        index=False, encoding='utf-8-sig')
    cases = pd.concat(rule_rows, ignore_index=True)
    cases.to_csv(OUT / 'V3_selected_CAS_rule_hypothetical_8h.csv',
                 index=False, encoding='utf-8-sig')

    summary = []
    for (temp, sc), d in frames.groupby(['temperature','scenario_name']):
        official = d[d['official_TWA_coverage_90pct']]
        summary.append({
            'temperature': temp, 'scenario_name': sc, 'all_documents': len(d),
            'official_OEL_90pct_eligible': len(official),
            'official_OEL_insufficient': len(d)-len(official),
            'hypo_8h_index_ge1_count': int(official['hypo_8h_TierA_index'].ge(1).sum()),
            'hypo_8h_index_ge10_count': int(official['hypo_8h_TierA_index'].ge(10).sum()),
            'hypo_8h_P5_ge1_count': int(official['hypo_8h_TierA_index_P5'].ge(1).sum()),
            'hypo_8h_P95_lt1_count': int(official['hypo_8h_TierA_index_P95'].lt(1).sum()),
            'index_ge1_unique_compositions': int(official.loc[
                official['hypo_8h_TierA_index'].ge(1), '조성지문'].nunique()),
            'ratio_is_illustrative': True
        })
    sums = pd.DataFrame(summary)
    sums.to_csv(OUT / 'V2_room_scenario_summary.csv', index=False,
                encoding='utf-8-sig')

    byrule = []
    for sc, d in cases.groupby('scenario_name'):
        d10 = d[d['hypo_8h_single_CAS_TierA_index'].ge(1)]
        byrule.append({
            'scenario_name': sc,
            'eligible_rule_CAS_rows': len(d),
            'single_CAS_hypo_8h_index_ge1': len(d10),
            'single_CAS_hypo_8h_ge1_documents': d10['문서키'].nunique(),
            'single_CAS_hypo_8h_ge10': int(d['hypo_8h_single_CAS_TierA_index'].ge(10).sum()),
            'interpretation': 'not real 8h measured exposure; source and room assumptions hypothetical'
        })
    pd.DataFrame(byrule).to_csv(
        OUT / 'V4_selected_CAS_room_scenario_summary.csv',
        index=False, encoding='utf-8-sig')

    meta = {
        'basis_documents': base['문서키'].nunique(),
        'temperatures': ['25C','35C'],
        'scenario_count': len(SETTINGS),
        'numerical_crosswalk_note': 'A threshold in product liquid mass % is NOT a TWA in air.',
        'main_caveat': 'Hypothetical finite-source 8h average model, NOT actual exposure or a compliance determination.',
        'parameters': scenario_inputs,
        'summary': summary,
        'CAS_specific_summary': byrule
    }
    (OUT / 'headline_ventilation_sensitivity.json').write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print('HYPOTHETICAL_ROOM_VENTILATION ' + json.dumps(meta, ensure_ascii=False))


if __name__ == '__main__':
    main()
