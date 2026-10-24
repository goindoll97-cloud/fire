# -*- coding: utf-8 -*-
"""Research-question tables and figures from v2/results/scenarios_270.csv.gz.

RQ1  Same-use variation of the mixture vapour hazard index (I_mix) and robust
     lower-hazard alternatives inside each use group.
RQ2  Effect of liquid-phase non-ideality: Raoult vs modified UNIFAC (Dortmund) vs HANNA.
RQ3  Robustness to MSDS composition ranges and UVCB surrogate composition.
RQ4  OEL coverage gap: share of saturated vapour without a Korean TWA, and how the
     gap changes product rankings (Tier A official only vs Tier B supplemented).
RQ5  Temperature sensitivity 25 -> 35 C.

Primary metric = HANNA, 25 C, Tier B index, central scenario; robustness from the
5th-95th percentile of sampled scenarios. All values are equilibrium head-space
potentials, not measured exposure.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

V2 = Path(__file__).resolve().parents[1]
RES = V2 / 'results'
USE_AUDIT = V2 / 'data' / 'use_classification_audit_270.csv'
OUT = RES / 'rq'
EXCLUDE_GROUPS = {'명칭만으로 판별 불가', '검색 오탐·비신나', ''}


def apply_use_audit(prod):
    """Apply documentary use review; fail closed on missing/outdated records.

    The original product classifications remain untouched in the 271-record source.
    An exploratory same-use comparison is NOT a validated field substitution.
    """
    if not USE_AUDIT.is_file():
        raise FileNotFoundError(f'270건 용도 감사표 누락: {USE_AUDIT}')
    aud = pd.read_csv(USE_AUDIT, dtype=str, encoding='utf-8-sig').fillna('')
    required = ['문서키', '기존제품ID', '원본SHA256', '용도군_감사', '비교세부군',
                '계절형_감사', '비교군_감사', '용도감사상태', '탐색비교허용', '현장대체검증']
    missing = [c for c in required if c not in aud]
    if missing:
        raise ValueError(f'용도 감사표 필수 열 누락: {missing}')
    if aud['문서키'].duplicated().any() or prod['문서키'].duplicated().any():
        raise ValueError('감사표 또는 모델 결과에 중복 문서키가 있습니다.')
    if set(aud['문서키']) != set(prod['문서키']):
        raise ValueError('감사표와 270건 모델 입력 문서키가 다릅니다. 용도감사를 갱신하세요.')
    if (aud['기존제품ID'].str.strip() == '').any():
        raise ValueError('용도 감사표에 제품ID 결측이 있습니다.')
    if not aud['탐색비교허용'].isin(['0', '1']).all():
        raise ValueError('탐색비교허용은 0 또는 1이어야 합니다.')
    if (aud['현장대체검증'] != '0').any():
        raise ValueError('TDS/공정 호환성 현장검증 전에 대체승인을 표시할 수 없습니다.')
    if '기존제품ID' in prod:
        chk = aud[['문서키', '기존제품ID']].merge(
            prod[['문서키', '기존제품ID']], on='문서키', suffixes=('_audit', '_model'), validate='one_to_one')
        if (chk['기존제품ID_audit'] != chk['기존제품ID_model']).any():
            raise ValueError('용도 감사표 제품 ID와 모델 입력 ID가 불일치합니다.')
    raw_path = V2.parent / 'data' / 'source_rebuild_20261023' / 'MSDS_원본대조_270.csv'
    if raw_path.is_file():
        raw = pd.read_csv(raw_path, dtype=str, encoding='utf-8-sig').fillna('')
        chk = aud[['문서키', '원본SHA256']].merge(
            raw[['문서키', '원본SHA256']], on='문서키', suffixes=('_audit', '_raw'), validate='one_to_one')
        if len(chk) != len(aud) or (chk['원본SHA256_audit'] != chk['원본SHA256_raw']).any():
            raise ValueError('감사표 SHA256과 원본 MSDS 문서 정보가 불일치합니다.')
    aud = aud[required].drop(columns=['기존제품ID', '원본SHA256'])
    out = prod.merge(aud, on='문서키', how='left', validate='one_to_one')
    out['용도군_원분류'] = out['용도군']
    out['용도군'] = out['용도군_감사']
    out['세부분류'] = out['비교세부군']
    out['탐색비교허용'] = out['탐색비교허용'].eq('1')
    out['현장대체검증'] = False
    out['비교군'] = out['비교군_감사']
    if out['용도군'].isna().any() or out['비교군'].isna().any():
        raise ValueError('감사표 병합 결과에 결측 분류가 있습니다.')
    return out

MIN_COVERAGE = 0.9   # a product is rankable only if >=90 % of its saturated vapour has an OEL (Tier B)
ALCOHOLS = {'67-56-1', '64-17-5', '67-63-0', '71-36-3', '78-83-1', '107-98-2', '111-76-2',
            '112-34-5', '123-42-2', '97-64-3', '2110-78-3'}


def q(s, p):
    return float(np.percentile(s, p))


def load():
    s = pd.read_csv(RES / 'scenarios_270.csv.gz')
    prod = pd.read_csv(RES / 'products_270.csv').fillna('')
    sp = pd.read_csv(RES / 'species_central_270.csv.gz')
    # composition-identical documents (same CAS and ranges) -> one representative
    prod['대표문서'] = prod.groupby('조성지문')['문서키'].transform('first')
    prod['중복여부'] = prod['문서키'] != prod['대표문서']
    # Evidence-checked same-use subgroup and season/grade are loaded from a
    # document-keyed audit, not inferred from the broad product family alone.
    prod = apply_use_audit(prod)
    return s, prod, sp


def summarise(s):
    keys = ['문서키', 'T', 'model']
    cen = s[s.scenario == 'central'].set_index(keys)
    smp = s[s.scenario != 'central'].groupby(keys)
    out = cen[['Ptot_kPa', 'Csat_total_ppm', 'I_mix_TierA', 'I_mix_TierB', 'OELcov_TierA',
               'OELcov_TierB', 'Hmax_CAS_TierA']].copy()
    for m in ('I_mix_TierA', 'I_mix_TierB'):
        out[m + '_P5'] = smp[m].quantile(0.05)
        out[m + '_P95'] = smp[m].quantile(0.95)
    return out.reset_index()


def rq1(summ, prod):
    base = summ[(summ['T'] == '25C') & (summ.model == 'HANNA')].merge(prod, on='문서키')
    base = base[(~base['중복여부']) & base['탐색비교허용'] & ~base['용도군'].isin(EXCLUDE_GROUPS)].copy()
    # same manufacturer + same product name = revisions/variants of one product, not alternatives
    base = base.sort_values('문서키').drop_duplicates(['제조사', '제품명'], keep='first')
    base['순위가능'] = base['OELcov_TierB'] >= MIN_COVERAGE
    rows, ranks = [], []
    for g, dall in base.groupby('비교군'):
        d = dall[dall['순위가능']].sort_values('I_mix_TierB')
        held = dall[~dall['순위가능']]
        if len(d) >= 2:
            lo, hi = d.iloc[0], d.iloc[-1]
            pairs = list(itertools.combinations(d.itertuples(index=False), 2))
            robust = sum(1 for a, b in pairs
                         if a.I_mix_TierB_P95 < b.I_mix_TierB_P5 or b.I_mix_TierB_P95 < a.I_mix_TierB_P5)
            rows.append({'비교군': g, '용도군': lo['용도군'], '순위가능 제품수': len(d), '판정유보(OEL공백) 제품수': len(held),
                         'I_min': round(lo.I_mix_TierB, 1), '최저제품': lo['제품명'],
                         'I_max': round(hi.I_mix_TierB, 1), '최고제품': hi['제품명'],
                         '최고/최저 배수': round(hi.I_mix_TierB / lo.I_mix_TierB, 2),
                         '중앙값': round(d.I_mix_TierB.median(), 1),
                         '제품쌍수': len(pairs), '조성불확실성하 분리쌍': robust,
                         '분리쌍비율': round(robust / len(pairs), 3) if pairs else np.nan})
        n = len(d)
        for r, row in enumerate(d.itertuples(index=False), 1):
            grade = '상대낮음' if r <= n / 3 else ('상대높음' if r > 2 * n / 3 else '상대중간')
            better = d[d.I_mix_TierB_P95 < row.I_mix_TierB_P5]
            alt = better.iloc[0] if len(better) else None
            ranks.append({'비교군': g, '순위(낮을수록 유리)': r, '제품': row.제품명, '제조사': row.제조사,
                          '기존제품ID': row.기존제품ID, 'I_mix_TierB(중앙)': round(row.I_mix_TierB, 1),
                          'P5': round(row.I_mix_TierB_P5, 1), 'P95': round(row.I_mix_TierB_P95, 1),
                          '상대등급': grade if n >= 3 else '비교군 소규모', 'OEL커버리지_A': round(row.OELcov_TierA, 3),
                          'OEL커버리지_B': round(row.OELcov_TierB, 3), '주요기여성분(TierA)': row.Hmax_CAS_TierA,
                          '확실히 더 낮은 대체제품': alt['제품명'] if alt is not None else '',
                          '대체시 저감배수': round(row.I_mix_TierB / alt.I_mix_TierB, 2) if alt is not None else ''})
        for row in held.itertuples(index=False):
            ranks.append({'비교군': g, '순위(낮을수록 유리)': '', '제품': row.제품명, '제조사': row.제조사,
                          '기존제품ID': row.기존제품ID, 'I_mix_TierB(중앙)': round(row.I_mix_TierB, 1),
                          'P5': round(row.I_mix_TierB_P5, 1), 'P95': round(row.I_mix_TierB_P95, 1),
                          '상대등급': '판정유보: 증기 중 노출기준 없는 성분 비율 높음',
                          'OEL커버리지_A': round(row.OELcov_TierA, 3), 'OEL커버리지_B': round(row.OELcov_TierB, 3),
                          '주요기여성분(TierA)': row.Hmax_CAS_TierA, '확실히 더 낮은 대체제품': '', '대체시 저감배수': ''})
    return pd.DataFrame(rows).sort_values('최고/최저 배수', ascending=False), pd.DataFrame(ranks)


def rq2(summ, prod, sp):
    c = summ[summ['T'] == '25C'].pivot_table(index='문서키', columns='model', values='I_mix_TierB')
    c = c.join(prod.set_index('문서키')[['제품명', '용도군', '비교군', '탐색비교허용', '중복여부', 'UNIFAC비고']])
    cen = sp[(sp['T'] == '25C') & (sp.model == 'HANNA')]
    alc = cen[cen.CAS.isin(ALCOHOLS)].groupby('문서키').x_liq.sum()
    c['알코올몰분율'] = alc.reindex(c.index).fillna(0)
    c['HANNA/Raoult'] = c['HANNA'] / c['Raoult']
    c['UNIFAC/Raoult'] = c['UNIFAC'] / c['Raoult']
    c['HANNA/UNIFAC'] = c['HANNA'] / c['UNIFAC']
    u = c[~c['중복여부']]
    summary = {
        'products': int(len(u)),
        'UNIFAC_computable': int(u['UNIFAC'].notna().sum()),
        'UNIFAC_missing_reason_counts': prod[~prod['중복여부']]['UNIFAC비고'].map(
            lambda s: s.split(':')[0] if s else 'OK').value_counts().to_dict(),
        'HANNA/Raoult_median': round(u['HANNA/Raoult'].median(), 3),
        'HANNA/Raoult_P5_P95': [round(q(u['HANNA/Raoult'], 5), 3), round(q(u['HANNA/Raoult'], 95), 3)],
        'HANNA/Raoult_max': round(u['HANNA/Raoult'].max(), 3),
        'share_HANNA/Raoult_gt_1.2': round(float((u['HANNA/Raoult'] > 1.2).mean()), 3),
        'HANNA/UNIFAC_median': round(u['HANNA/UNIFAC'].median(), 3),
        'HANNA/UNIFAC_P5_P95': [round(q(u['HANNA/UNIFAC'].dropna(), 5), 3),
                                round(q(u['HANNA/UNIFAC'].dropna(), 95), 3)],
        'spearman_HANNA/Raoult_vs_alcohol_x': round(u[['HANNA/Raoult', '알코올몰분율']].corr('spearman').iloc[0, 1], 3),
    }
    rev = []
    g = u[u['탐색비교허용'] & ~u['용도군'].isin(EXCLUDE_GROUPS)]
    for grp, d in g.groupby('비교군'):
        for (ia, a), (ib, b) in itertools.combinations(d.iterrows(), 2):
            if np.sign(a.Raoult - b.Raoult) != np.sign(a.HANNA - b.HANNA):
                rev.append({'용도군': grp, '제품A': a.제품명, '제품B': b.제품명,
                            'Raoult A/B': round(a.Raoult / b.Raoult, 3), 'HANNA A/B': round(a.HANNA / b.HANNA, 3),
                            '알코올몰분율 A': round(a.알코올몰분율, 3), '알코올몰분율 B': round(b.알코올몰분율, 3)})
    npairs = sum(len(d) * (len(d) - 1) // 2 for _, d in g.groupby('비교군'))
    summary['within_group_pairs'] = int(npairs)
    summary['raoult_vs_hanna_order_reversals'] = len(rev)
    return c.reset_index(), pd.DataFrame(rev), summary


def rq3(summ, prod):
    b = summ[(summ['T'] == '25C') & (summ.model == 'HANNA')].merge(prod, on='문서키')
    b = b[~b['중복여부']]
    b['상대폭(P95/P5)'] = b.I_mix_TierB_P95 / b.I_mix_TierB_P5
    return b[['문서키', '제품명', '용도군', 'I_mix_TierB', 'I_mix_TierB_P5', 'I_mix_TierB_P95', '상대폭(P95/P5)',
              'UVCB대리물질사용', '조성비고']], {
        'relative_width_median': round(b['상대폭(P95/P5)'].median(), 3),
        'relative_width_P90': round(q(b['상대폭(P95/P5)'], 90), 3)}


def rq4(summ, prod, sp):
    b = summ[(summ['T'] == '25C') & (summ.model == 'HANNA')].merge(prod, on='문서키')
    b = b[~b['중복여부']]
    cen = sp[(sp['T'] == '25C') & (sp.model == 'HANNA')]
    oel = pd.read_csv(V2 / 'data' / 'oel_assignment.csv')
    noA = set(oel[oel.TierA_TWA_ppm.isna()].CAS)
    gap = cen[cen.CAS.isin(noA) & cen.문서키.isin(b.문서키)]
    contrib = gap.groupby(['CAS', 'name']).agg(제품수=('문서키', 'nunique'), 증기몰분율합=('y_vap', 'sum'),
                                               증기몰분율최대=('y_vap', 'max')).reset_index()
    contrib = contrib.sort_values('증기몰분율합', ascending=False)
    # rank shift within groups: Tier A (official only) vs Tier B
    shifts = []
    for g, d in b[b['탐색비교허용'] & ~b['용도군'].isin(EXCLUDE_GROUPS)].groupby('비교군'):
        if len(d) < 3:
            continue
        ra, rb = d.I_mix_TierA.rank(), d.I_mix_TierB.rank()
        for (_, row), a_, b_ in zip(d.iterrows(), ra, rb):
            shifts.append({'용도군': g, '제품': row.제품명, 'TierA순위': int(a_), 'TierB순위': int(b_),
                           '순위하락(공식기준에서 과소평가)': int(b_ - a_), 'OEL커버리지_A': round(row.OELcov_TierA, 3)})
    shifts = pd.DataFrame(shifts)
    summary = {
        'products': int(len(b)),
        'coverageA_median': round(b.OELcov_TierA.median(), 3),
        'n_coverageA_lt_0.9': int((b.OELcov_TierA < 0.9).sum()),
        'n_coverageA_lt_0.5': int((b.OELcov_TierA < 0.5).sum()),
        'coverageB_median': round(b.OELcov_TierB.median(), 3),
        'n_coverageB_lt_0.9': int((b.OELcov_TierB < 0.9).sum()),
        'n_products_rank_worse_under_B_by_ge2': int((shifts['순위하락(공식기준에서 과소평가)'] >= 2).sum()) if len(shifts) else 0,
        'IB_over_IA_median': round(float((b.I_mix_TierB / b.I_mix_TierA).median()), 3),
        'IB_over_IA_max': round(float((b.I_mix_TierB / b.I_mix_TierA).replace(np.inf, np.nan).max()), 3),
    }
    return b[['문서키', '제품명', '용도군', 'OELcov_TierA', 'OELcov_TierB', 'I_mix_TierA', 'I_mix_TierB']], contrib, shifts, summary


def rq5(summ, prod):
    p = summ[summ.model == 'HANNA'].pivot_table(index='문서키', columns='T', values=['I_mix_TierB', 'Csat_total_ppm'])
    r = (p['I_mix_TierB']['35C'] / p['I_mix_TierB']['25C']).rename('I(35C)/I(25C)')
    r = r.to_frame().join(prod.set_index('문서키')[['제품명', '용도군', '계절변형', '중복여부']])
    u = r[~r['중복여부']]
    return u.reset_index(), {'ratio_median': round(u['I(35C)/I(25C)'].median(), 3),
                             'ratio_range': [round(u['I(35C)/I(25C)'].min(), 3), round(u['I(35C)/I(25C)'].max(), 3)]}


def figures(summ, prod, rq1_rank, rq2_tab):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    plt.rcParams['font.family'] = ['DejaVu Sans']
    # Fig 1: within-group distribution (log) of I_mix Tier B, HANNA 25 C
    d = rq1_rank.copy()
    d = d[d['순위(낮을수록 유리)'] != '']
    d = d[d.groupby('비교군')['제품'].transform('count') >= 2]
    groups = d.groupby('비교군')['I_mix_TierB(중앙)'].median().sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for i, g in enumerate(groups):
        v = d[d['비교군'] == g]
        ax.errorbar(np.full(len(v), i) + np.linspace(-0.2, 0.2, len(v)), v['I_mix_TierB(중앙)'],
                    yerr=[v['I_mix_TierB(중앙)'] - v['P5'], v['P95'] - v['I_mix_TierB(중앙)']],
                    fmt='o', ms=3, lw=0.6, color='#2a6f97', alpha=0.85)
    ax.set_yscale('log')
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f'G{i + 1}' for i in range(len(groups))], rotation=90, fontsize=8)
    ax.set_ylabel('Mixture vapour hazard index (Tier B, HANNA, 25 C)')
    ax.set_xlabel('Comparison group: use x drying/season grade (legend in fig1_groups.csv)')
    fig.tight_layout()
    fig.savefig(OUT / 'fig1_within_group_index.png', dpi=200)
    pd.DataFrame({'code': [f'G{i + 1}' for i in range(len(groups))], '비교군': groups}).to_csv(
        OUT / 'fig1_groups.csv', index=False, encoding='utf-8-sig')
    plt.close(fig)
    # Fig 2: HANNA/Raoult vs alcohol mole fraction
    t = rq2_tab[~rq2_tab['중복여부']]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(t['알코올몰분율'], t['HANNA/Raoult'], s=10, color='#c1121f', label='HANNA / Raoult')
    ax.scatter(t['알코올몰분율'], t['UNIFAC/Raoult'], s=10, marker='x', color='#555', label='UNIFAC-Do / Raoult')
    ax.axhline(1, lw=0.8, color='k')
    ax.set_xlabel('Alcohol / glycol-ether mole fraction in liquid')
    ax.set_ylabel('Index ratio vs ideal solution')
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'fig2_nonideality.png', dpi=200)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s, prod, sp = load()
    summ = summarise(s)
    summ.to_csv(OUT / 'product_summary_all_models.csv', index=False, encoding='utf-8-sig')
    r1, r1rank = rq1(summ, prod)
    r1.to_csv(OUT / 'RQ1_group_fold_range.csv', index=False, encoding='utf-8-sig')
    r1rank.to_csv(OUT / 'RQ1_field_ranking_table.csv', index=False, encoding='utf-8-sig')
    r2tab, r2rev, r2sum = rq2(summ, prod, sp)
    r2tab.to_csv(OUT / 'RQ2_model_ratios.csv', index=False, encoding='utf-8-sig')
    r2rev.to_csv(OUT / 'RQ2_rank_reversals.csv', index=False, encoding='utf-8-sig')
    r3tab, r3sum = rq3(summ, prod)
    r3tab.to_csv(OUT / 'RQ3_uncertainty.csv', index=False, encoding='utf-8-sig')
    r4tab, r4con, r4shift, r4sum = rq4(summ, prod, sp)
    r4tab.to_csv(OUT / 'RQ4_oel_coverage.csv', index=False, encoding='utf-8-sig')
    r4con.to_csv(OUT / 'RQ4_uncovered_species.csv', index=False, encoding='utf-8-sig')
    r4shift.to_csv(OUT / 'RQ4_rank_shift_tierA_vs_B.csv', index=False, encoding='utf-8-sig')
    r5tab, r5sum = rq5(summ, prod)
    r5tab.to_csv(OUT / 'RQ5_temperature.csv', index=False, encoding='utf-8-sig')
    figures(summ, prod, r1rank, r2tab)
    headline = {
        'population': {'documents': int(len(prod)), 'composition_unique': int((~prod['중복여부']).sum()),
                       'with_UVCB_surrogates': int(prod['UVCB대리물질사용'].astype(str).eq('True').sum()),
                       'range_repaired': int((prod['조성비고'] != '').sum())},
        'RQ1': {'comparison_groups_with_2plus': int(len(r1)),
                'rankable_products': int(r1['순위가능 제품수'].sum()),
                'held_for_OEL_gap': int((r1rank['상대등급'].astype(str).str.startswith('판정유보')).sum()),
                'fold_range_median': round(float(r1['최고/최저 배수'].median()), 2),
                'fold_range_IQR': [round(float(r1['최고/최저 배수'].quantile(.25)), 2), round(float(r1['최고/최저 배수'].quantile(.75)), 2)],
                'fold_range_max': round(float(r1['최고/최저 배수'].max()), 2),
                'robust_pair_share_overall': round(float(r1['조성불확실성하 분리쌍'].sum() / r1['제품쌍수'].sum()), 3),
                'products_with_robust_lower_alternative': int((r1rank['확실히 더 낮은 대체제품'].astype(str) != '').sum()),
                'median_reduction_factor_if_switched': round(float(pd.to_numeric(r1rank['대체시 저감배수'], errors='coerce').median()), 2)},
        'RQ2': r2sum, 'RQ3': r3sum, 'RQ4': r4sum, 'RQ5': r5sum}
    (OUT / 'headline_numbers.json').write_text(json.dumps(headline, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(headline, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
