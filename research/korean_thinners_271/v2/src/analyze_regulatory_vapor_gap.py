# -*- coding: utf-8 -*-
"""Preliminary CAS-specific regulatory *mass fraction* vs saturated-vapor *OEL* audit.

WARNING:
- The official-law list is NOT comprehensive: selected KREACH-verified CAS rules only.
- Below threshold of one provision DOES NOT equal unregulated or safe under any law.
- Product mass percentage is legally distinct from saturated vapor concentration ppm.
- Vapor calculations (25C HANNA) are equilibrium potential, NOT workplace 8h TWA.
- The original Section 3 composition ranges are still unverified extracts.
- Xylene vapor from surrogate naphtha cannot always be separated from declared xylene.
- Hazardous-substance law effective-date exceptions must be individually checked.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import analyze_v2 as a2

V2 = a2.V2
BASE = V2.parent
OUT = a2.RES / 'regulatory_gap'
COMPOSITION = BASE / 'data' / 'source_rebuild_20261023' / 'MSDS_구성성분_1454.csv'
RULES = V2 / 'data' / 'verified_regulatory_thresholds_20261024.csv'
SPECIES = a2.RES / 'species_central_270.csv.gz'
P_ATM = 101.325

# These are actual model species representing the ORIGINAL CAS, not
# a direct legal classification of vapor mixture.
SPECIES_MAP = {
    '1330-20-7': ('108-38-3', '106-42-3', '95-47-6'),
}
# Xylene-model isomers can ALSO come from an aromatic-naphtha UVCB component.
CROSS_ORIGIN = {'1330-20-7': {'64742-95-6'}}


def status(lo, hi, threshold):
    if not np.isfinite(lo) or not np.isfinite(hi) or lo < 0 or hi > 100 or hi < lo:
        return 'range_invalid'
    if hi < threshold:
        return 'below_rule_threshold'
    if lo >= threshold:
        return 'at_or_above_rule_threshold'
    return 'straddles_threshold'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    comp = pd.read_csv(COMPOSITION, encoding='utf-8-sig', dtype={'CAS': str})
    rules = pd.read_csv(RULES, encoding='utf-8-sig', dtype={'CAS': str})
    _, products, _ = a2.load()
    species = pd.read_csv(SPECIES, encoding='utf-8-sig', dtype={'CAS': str})
    oel = pd.read_csv(V2 / 'data' / 'oel_assignment.csv', encoding='utf-8-sig', dtype={'CAS': str})
    properties = pd.read_csv(V2 / 'data' / 'pure_properties.csv', encoding='utf-8-sig', dtype={'CAS': str})
    mw = properties.set_index('CAS')['MW_g_mol'].astype(float).to_dict()
    surr = pd.read_csv(V2 / 'data' / 'surrogates.csv', encoding='utf-8-sig', dtype=str)
    surrogate_components = surr.groupby('UVCB_CAS')['component_CAS'].agg(set).to_dict()
    oel_a = oel.set_index('CAS')['TierA_TWA_ppm'].to_dict()

    main_sp = species[(species['T'] == '25C')].copy()
    docs = products[['문서키', '기존제품ID', '제조사', '제품명', '조성지문',
                     'UVCB대리물질사용', '용도군', '탐색비교허용']].copy()
    cs = comp[['문서키', 'CAS', '함량하한_pct', '함량상한_pct']].copy()
    cs['CAS'] = cs['CAS'].astype(str).str.strip()
    # CAS can occur in multiple rows of a single source MSDS: aggregate the
    # bounds before threshold comparison. Groupby min/max is not sufficient.
    cs = cs.groupby(['문서키', 'CAS'], as_index=False)[
        ['함량하한_pct', '함량상한_pct']].sum(min_count=1)
    declared = cs.merge(rules, on='CAS', how='inner', validate='many_to_many')
    declared = declared.merge(docs, on='문서키', how='left', validate='many_to_one')
    if declared['기존제품ID'].isna().any():
        raise ValueError('Source constituent row could not be matched to 270 documents.')
    declare_cas = comp.groupby('문서키')['CAS'].agg(set).to_dict()
    sp_groups = dict(tuple(main_sp.groupby(['문서키', 'model'])))
    index_counts = {}
    output = []
    for r in declared.itertuples(index=False):
        cas = str(r.CAS)
        rule_name = str(r.regulation_category)
        threshold = float(r.threshold_mass_pct)
        lo = float(r.함량하한_pct)
        hi = float(r.함량상한_pct)
        legal = status(lo, hi, threshold)
        model_sp = SPECIES_MAP.get(cas, (cas,))
        # Overlap from a listed UVCB is a provenance uncertainty; do NOT
        # credit ALL species vapor to declared CAS if it also came from UVCB.
        origin_overlap = any(
            cas != uv and any(s in set(surrogate_components.get(uv, set())) for s in model_sp)
            for uv in declare_cas.get(r.문서키, set())
        )
        for model in ('HANNA', 'Raoult', 'UNIFAC'):
            d = sp_groups.get((r.문서키, model))
            sub = d[d['CAS'].isin(model_sp)] if d is not None else pd.DataFrame()
            if sub.empty:
                vapor_ppm = np.nan
                term_a = np.nan
                vapor_y = np.nan
                liquid_x = np.nan
                ps = np.nan
                vwt = np.nan
                check = 'model_or_species_not_computed'
            else:
                ps = float(sub['p_kPa'].sum())
                vapor_ppm = ps / P_ATM * 1e6
                term_a = float(sub['term_TierA'].sum())
                vapor_y = float(sub['y_vap'].sum())
                liquid_x = float(sub['x_liq'].sum())
                check = 'surrogate_origin_overlap' if origin_overlap else 'central_model_value'
                # Mass fraction in the solvent vapor (air excluded).
                # Include all available species for correct vapor denominator.
                allweights = d['p_kPa'] * d['CAS'].map(mw)
                numer = float((sub['p_kPa'] * sub['CAS'].map(mw)).sum())
                vwt = 100 * numer / float(allweights.sum()) if float(allweights.sum()) > 0 else np.nan
            index_counts[(cas, model)] = index_counts.get((cas, model), 0) + 1
            output.append({
                '문서키': r.문서키, '기존제품ID': r.기존제품ID,
                '제조사': r.제조사, '제품명': r.제품명, '용도군': r.용도군,
                '비교용도감사허용': r.탐색비교허용,
                '조성지문': r.조성지문, 'UVCB대리물질사용': r.UVCB대리물질사용,
                'CAS_원문': cas, '물질명': r.name_en,
                '규정범주': rule_name, '고유번호': r.item_number,
                '규정_액상함량기준_mass_pct': threshold,
                'MSDS_액상하한_mass_pct': lo, 'MSDS_액상상한_mass_pct': hi,
                '해당규정_함량비교': legal, '기상모델': model,
                '모형CAS': '|'.join(model_sp),
                '기상성분추적신뢰': check,
                '액상분자몰분율_x_중앙조성': liquid_x,
                '기상용제내_몰분율_y': vapor_y,
                '기상용제내_질량분율_pct': vwt,
                '모형_기상분압_kPa': ps,
                '모형_포화기상농도_ppm_in_air': vapor_ppm,
                '모형_개별성분_TierA_OEL지수': term_a,
                'OEL_공식CAS근거': '|'.join(
                    f'{x}={oel_a.get(x, np.nan)}' for x in model_sp),
                '규정출처': r.source_url,
                '법적취급': '이 단일 규정의 혼합물 액상함량 기준만 비교. 기타 법령 적용 여부 별도 확인',
                '포화모형주의': '포화평형 잠재농도이며 실제 개인 노출농도나 TWA초과 아님'
            })
    df = pd.DataFrame(output)
    if len(df) == 0:
        raise ValueError('No selected CAS regulatory clauses matched the source MSDS data')
    df.to_csv(OUT / 'R1_constituent_rule_vapor_crosswalk.csv',
              index=False, encoding='utf-8-sig')
    primary = df[df['기상모델'] == 'HANNA'].copy()
    primary['함량미만_성분포화지수10'] = (
        primary['해당규정_함량비교'].eq('below_rule_threshold') &
        primary['모형_개별성분_TierA_OEL지수'].ge(10) &
        primary['기상성분추적신뢰'].eq('central_model_value'))
    primary['함량미만_성분포화지수100'] = (
        primary['해당규정_함량비교'].eq('below_rule_threshold') &
        primary['모형_개별성분_TierA_OEL지수'].ge(100) &
        primary['기상성분추적신뢰'].eq('central_model_value'))
    primary['원문CAS별_함량아래_포화10_추적확인'] = primary['함량미만_성분포화지수10']
    grouped = primary.groupby(['CAS_원문', '물질명', '규정범주', '규정_액상함량기준_mass_pct'],
                              dropna=False).agg(
        문서수=('문서키', 'nunique'),
        기준미만문서수=('해당규정_함량비교', lambda s: int((s == 'below_rule_threshold').sum())),
        기준이상문서수=('해당규정_함량비교', lambda s: int((s == 'at_or_above_rule_threshold').sum())),
        기준경계문서수=('해당규정_함량비교', lambda s: int((s == 'straddles_threshold').sum())),
        출처겹침문서수=('기상성분추적신뢰', lambda s: int((s == 'surrogate_origin_overlap').sum())),
        기준미만_포화10이상_추적확인=('함량미만_성분포화지수10', 'sum'),
        기준미만_포화100이상_추적확인=('함량미만_성분포화지수100', 'sum')
    ).reset_index()
    grouped.to_csv(OUT / 'R2_regulatory_rule_summary.csv', index=False,
                   encoding='utf-8-sig')

    # Scientific phase-enrichment test: both sides are mass percent, but the
    # *liquid product* legal cutoff MUST NOT be applied to a *solvent-only vapor*.
    # This is a hypothetical numeric boundary crossing, not legal exceedance.
    df['기상질량분율_가상경계교차_법적해당아님'] = (
        df['해당규정_함량비교'].eq('below_rule_threshold') &
        df['기상용제내_질량분율_pct'].ge(df['규정_액상함량기준_mass_pct']) &
        df['기상성분추적신뢰'].eq('central_model_value'))
    df['액상대비기상_몰분율농축배율'] = np.where(
        df['액상분자몰분율_x_중앙조성'] > 0,
        df['기상용제내_몰분율_y'] / df['액상분자몰분율_x_중앙조성'],
        np.nan)
    # Initial R1 includes source-trace audit, now with phase-specific flags.
    df.to_csv(OUT / 'R1_constituent_rule_vapor_crosswalk.csv',
              index=False, encoding='utf-8-sig')
    vapor_main = df[df['기상모델'].eq('HANNA')].copy()
    vapor_main['가상경계교차_예비'] = vapor_main['기상질량분율_가상경계교차_법적해당아님']
    vgroups = vapor_main.groupby(
        ['CAS_원문', '물질명', '규정범주',
         '규정_액상함량기준_mass_pct'], dropna=False).agg(
            액상기준미만문서수=('해당규정_함량비교',
                             lambda s: int(s.eq('below_rule_threshold').sum())),
            기상_가상질량분율경계교차_법적비해당=('가상경계교차_예비', 'sum'),
            기상성분추적불확실_문서수=('기상성분추적신뢰',
                                  lambda s: int(s.eq('surrogate_origin_overlap').sum())),
            액상기상_몰분율농축배율_중앙값=('액상대비기상_몰분율농축배율', 'median')
        ).reset_index()
    vgroups.to_csv(OUT / 'R6_vapor_only_phase_enrichment.csv',
                   index=False, encoding='utf-8-sig')
    phase_cases = vapor_main[vapor_main['가상경계교차_예비']].copy()
    phase_cases[['문서키', '기존제품ID', '제품명', 'CAS_원문', '규정범주',
                 'MSDS_액상하한_mass_pct', 'MSDS_액상상한_mass_pct',
                 '규정_액상함량기준_mass_pct', '기상용제내_질량분율_pct',
                 '액상대비기상_몰분율농축배율',
                 '모형_개별성분_TierA_OEL지수']].to_csv(
        OUT / 'R7_hypothetical_phase_fraction_cases_not_legal.csv',
        index=False, encoding='utf-8-sig')


    # Three activity-coefficient baselines are compared on the SAME source
    # CAS / document / regulation clause, not on different model populations.
    mod_keys = ['문서키', 'CAS_원문', '규정범주', '고유번호']
    cv = df.pivot_table(index=mod_keys, columns='기상모델',
                        values='모형_개별성분_TierA_OEL지수', aggfunc='first')
    elig = primary[primary['해당규정_함량비교'].eq('below_rule_threshold') &
                   primary['기상성분추적신뢰'].eq('central_model_value')]
    cv = cv.join(elig.set_index(mod_keys)[
        ['해당규정_함량비교']], how='inner')
    common_models = cv.dropna(subset=['HANNA', 'Raoult', 'UNIFAC']).copy()
    for model in ('HANNA', 'Raoult', 'UNIFAC'):
        common_models[model + '_over_TWA_10x'] = common_models[model].ge(10)
    common_models.to_csv(OUT / 'R4_common_case_model_sensitivity.csv',
                         encoding='utf-8-sig')
    model_robustness = {
        'eligible_common_three_model_cases': int(len(common_models)),
        'all_three_over_10x': int(common_models[[
            'HANNA_over_TWA_10x', 'Raoult_over_TWA_10x',
            'UNIFAC_over_TWA_10x']].all(axis=1).sum()),
        'any_three_models_over_10x': int(common_models[[
            'HANNA_over_TWA_10x', 'Raoult_over_TWA_10x',
            'UNIFAC_over_TWA_10x']].any(axis=1).sum()),
    }
    matched = primary[primary['함량미만_성분포화지수10']].copy()
    examples = matched.sort_values('모형_개별성분_TierA_OEL지수', ascending=False)
    examples = examples.drop_duplicates(['CAS_원문', '고유번호']).head(8)
    examples = examples[[
        '기존제품ID', '제품명', 'CAS_원문', '규정범주',
        'MSDS_액상하한_mass_pct', 'MSDS_액상상한_mass_pct',
        '규정_액상함량기준_mass_pct', '모형_포화기상농도_ppm_in_air',
        '모형_개별성분_TierA_OEL지수', 'UVCB대리물질사용']]
    bydoc = primary.groupby('문서키').agg(
        검토규정_연결성분수=('CAS_원문', 'nunique'),
        비교된규정행수=('CAS_원문', 'size'),
        일부규정_함량미만포화10이상_존재=('함량미만_성분포화지수10', 'max'),
        일부규정_함량미만포화100이상_존재=('함량미만_성분포화지수100', 'max')
    ).reset_index()
    bydoc = products[['문서키', '기존제품ID', '제조사', '제품명',
                      '조성지문', 'UVCB대리물질사용']].merge(bydoc, on='문서키', how='left')
    for c in ('일부규정_함량미만포화10이상_존재', '일부규정_함량미만포화100이상_존재'):
        bydoc[c] = bydoc[c].fillna(False).astype(bool)
    bydoc.to_csv(OUT / 'R3_270_document_coverage.csv', index=False,
                 encoding='utf-8-sig')

    # Isolate samples of 270 docs, NOT iid exposure measurements.
    counts = {
        'documents_total': int(len(bydoc)),
        'documents_with_any_verified_CAS_rule': int(bydoc['검토규정_연결성분수'].notna().sum()),
        'cas_rule_rows_HANNA': int(len(primary)),
        'unique_CAS_reviewed': sorted(primary['CAS_원문'].unique().tolist()),
        'source_CAS_total': int(comp['CAS'].nunique()),
        'CAS_with_reviewed_rules': int(primary['CAS_원문'].nunique()),
        'cases_below_clause_and_single_TierA_index_ge_10_traceable':
            int(primary['함량미만_성분포화지수10'].sum()),
        'documents_with_at_least_one_case_ge_10':
            int(bydoc['일부규정_함량미만포화10이상_존재'].sum()),
        'cases_below_clause_and_single_TierA_index_ge_100_traceable':
            int(primary['함량미만_성분포화지수100'].sum()),
        'documents_with_at_least_one_case_ge_100':
            int(bydoc['일부규정_함량미만포화100이상_존재'].sum()),
        'rule_status_distribution': primary['해당규정_함량비교'].value_counts().to_dict(),
        'composition_unique_in_total': int(products['조성지문'].nunique()),
        'documents_flagged_unique_composition': int(
            bydoc.loc[bydoc['일부규정_함량미만포화10이상_존재'],
                      '조성지문'].nunique()),
        'vapor_only_crossing_documents_unique_not_legal':
            int(phase_cases['문서키'].nunique()),
        'vapor_only_crossing_compositions_unique_not_legal':
            int(phase_cases['조성지문'].nunique()),
        'vapor_only_mass_pct_hypothetical_crossing_count_not_legal':
            int(vapor_main['가상경계교차_예비'].sum()),
        'vapor_only_numeric_boundary_crossings_by_CAS': vgroups[
            ['CAS_원문', '규정범주',
             '기상_가상질량분율경계교차_법적비해당']].to_dict(orient='records'),
        'model_crosscheck': model_robustness,
        'with_UVCB_origin_overlap_cases': int(
            primary['기상성분추적신뢰'].eq('surrogate_origin_overlap').sum()),
        'note_1': 'Excluded below-threshold flag does not determine other-acts legal applicability.',
        'note_2': 'The numerator is saturated equilibrium ppm/OEL, not workplace exposure or TWA exceedance.',
        'note_3': 'This is a SCREENING of selected verified CAS rules, not a complete 2026 designation annex.',
        'note_4': 'No benzene was directly listed in the 270-MSDS Section 3 extracted CAS rows.',
        'note_5': 'Xylene isomer vapor may partly derive from aromatic UVCB proxy rather than the declared xylene alone.',
    }
    (OUT / 'headline_regulatory_vapor_gap.json').write_text(
        json.dumps(counts, ensure_ascii=False, indent=2), encoding='utf-8')
    phase_samples = phase_cases.sort_values('기상용제내_질량분율_pct', ascending=False)
    phase_samples = phase_samples.drop_duplicates(['CAS_원문', '규정범주']).head(7)
    print('REGULATORY_VAPOR_PHASE_CASES ' +
          json.dumps(phase_samples[[
              '기존제품ID', '제품명', 'CAS_원문', '규정범주',
              'MSDS_액상하한_mass_pct', 'MSDS_액상상한_mass_pct',
              '규정_액상함량기준_mass_pct', '기상용제내_질량분율_pct',
              '액상대비기상_몰분율농축배율']].to_dict(orient='records'), ensure_ascii=False))
    print('REGULATORY_VAPOR_ENRICHMENT ' +
          json.dumps(vgroups.to_dict(orient='records'), ensure_ascii=False))
    print('REGULATORY_VAPOR_RULE_GROUPS ' +
          json.dumps(grouped.to_dict(orient='records'), ensure_ascii=False))
    print('REGULATORY_VAPOR_EXAMPLES ' +
          json.dumps(examples.to_dict(orient='records'), ensure_ascii=False))
    print('REGULATORY_VAPOR_GAP ' + json.dumps(counts, ensure_ascii=False))


if __name__ == '__main__':
    main()
