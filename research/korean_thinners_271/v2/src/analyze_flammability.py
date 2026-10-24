# -*- coding: utf-8 -*-
"""RQ6–RQ8: 화재 지표 검증과 건강·화재 동시 평가.

RQ6  예측 인화점(Liaw)을 MSDS 제9항 인화점과 비교한다 (Raoult / UNIFAC / HANNA).
RQ7  같은 용도·계절군 안에서 흡입 노출 지수(I_mix)와 인화점이 함께 움직이는가?
     (Spearman, 사분면, 두 위험이 엇갈리는 제품)
RQ8  같은 비교군 안에서 노출 지수는 확실히 낮고 인화점은 확실히 높은 대체 제품이 있는가?
     (위험물안전관리법 석유류 구분이 바뀌는 경우 포함)

주지표: HANNA, 중앙조성, 25°C, Tier B. 불확실성은 샘플 시나리오의 P5–P95.
인화점 구분 기준(위험물안전관리법 시행령 별표 1 제4류): 제1석유류 21°C 미만, 제2석유류 21°C 이상 70°C 미만.
지정수량: 비수용성 제1석유류 200 L, 제2석유류 1,000 L. 투고 전 법제처 현행 조문으로 확인한다.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

import analyze_v2 as a2

V2 = a2.V2
RES = a2.RES
OUT = RES / 'rq'
CLASS1, CLASS2 = 21.0, 70.0


def petroleum_class(t):
    if pd.isna(t):
        return ''
    return '제1석유류' if t < CLASS1 else ('제2석유류' if t < CLASS2 else '제3석유류 이상')


def fire_summary(f):
    keys = ['문서키', 'model']
    cen = f[f.scenario == 'central'].set_index(keys)
    smp = f[f.scenario != 'central'].groupby(keys)
    out = cen[['Tfp_Liaw_C', 'F25', 'F35', 'LELmix25_pct', 'Csat25_pct', 'FP추정성분포함']].copy()
    out['Tfp_P5'] = smp['Tfp_Liaw_C'].quantile(0.05)
    out['Tfp_P95'] = smp['Tfp_Liaw_C'].quantile(0.95)
    return out.reset_index()


def rq6(fs, prod, msds):
    d = fs.merge(prod[['문서키', '기존제품ID', '중복여부', '용도군']], on='문서키').merge(
        msds[['기존제품ID', '인화점_C']].rename(columns={'인화점_C': 'FP_MSDS'}), on='기존제품ID')
    d = d[~d['중복여부']].dropna(subset=['Tfp_Liaw_C', 'FP_MSDS'])
    rows = []
    for m, g in d.groupby('model'):
        err = g['Tfp_Liaw_C'] - g['FP_MSDS']
        inside = ((g['FP_MSDS'] >= g['Tfp_P5']) & (g['FP_MSDS'] <= g['Tfp_P95'])).mean()
        rows.append({'model': m, 'n': len(g), 'MAE_C': round(err.abs().mean(), 2), 'median_err_C': round(err.median(), 2),
                     'bias_C': round(err.mean(), 2), 'within_5C': round((err.abs() <= 5).mean(), 3),
                     'within_10C': round((err.abs() <= 10).mean(), 3),
                     'spearman': round(g['Tfp_Liaw_C'].corr(g['FP_MSDS'], method='spearman'), 3),
                     'MSDS값이 시나리오 P5–P95 안': round(inside, 3),
                     '21C 구분 일치율': round(((g['Tfp_Liaw_C'] < CLASS1) == (g['FP_MSDS'] < CLASS1)).mean(), 3),
                     '23C 구분 일치율': round(((g['Tfp_Liaw_C'] < 23) == (g['FP_MSDS'] < 23)).mean(), 3)})
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / 'RQ6_flashpoint_validation.csv', index=False, encoding='utf-8-sig')
    h = d[d.model == 'HANNA'].copy()
    h['err'] = h['Tfp_Liaw_C'] - h['FP_MSDS']
    h.sort_values('err').to_csv(OUT / 'RQ6_flashpoint_by_product.csv', index=False, encoding='utf-8-sig')
    return tab, h


def rq6_bounds(prod, msds, sp):
    """모델과 무관한 일관성 검사: 이상용액에서 혼합물 인화점은 성분 인화점의 최솟값과 최댓값 사이에 있다.
    MSDS 인화점이 이 범위를 벗어나면 조성과 모순된다. 또 몰분율 5% 이상인 가장 낮은 인화점 성분보다 25°C 넘게 높은 값도 센다."""
    fl = pd.read_csv(V2 / 'data' / 'pure_flammability.csv').set_index('CAS')['FP_C']
    s = sp[(sp.model == 'Raoult') & (sp['T'] == '25C')].copy()
    s['fp'] = s['CAS'].map(fl)
    lo = s.groupby('문서키')['fp'].min().rename('FP_min')
    hi = s.groupby('문서키')['fp'].max().rename('FP_max')
    lo5 = s[s.x_liq >= 0.05].groupby('문서키')['fp'].min().rename('FP_min_5mol')
    d = prod[~prod['중복여부']][['문서키', '기존제품ID', '제품명']].merge(
        msds[['기존제품ID', '인화점_C']].rename(columns={'인화점_C': 'FP_MSDS'}), on='기존제품ID')
    d = d.join(lo, on='문서키').join(hi, on='문서키').join(lo5, on='문서키')
    d['below_min'] = d['FP_MSDS'] < d['FP_min']
    d['above_max'] = d['FP_MSDS'] > d['FP_max']
    d['over_lowest_5mol_by_25C'] = d['FP_MSDS'] > d['FP_min_5mol'] + 25
    d.to_csv(OUT / 'RQ6_msds_flash_consistency.csv', index=False, encoding='utf-8-sig')
    return {'consistency_n': int(len(d)), 'msds_below_lowest_component_fp': int(d['below_min'].sum()),
            'msds_above_highest_component_fp': int(d['above_max'].sum()),
            'msds_over_lowest_5mol_component_by_25C': int(d['over_lowest_5mol_by_25C'].sum()),
            'msds_top_values': msds['인화점_C'].value_counts().head(5).to_dict()}


def rq7_rq8(fs, summ, prod):
    health = summ[(summ['T'] == '25C') & (summ.model == 'HANNA')][
        ['문서키', 'I_mix_TierB', 'I_mix_TierB_P5', 'I_mix_TierB_P95', 'OELcov_TierB']]
    fire = fs[fs.model == 'HANNA'].drop(columns='model')
    b = prod.merge(health, on='문서키').merge(fire, on='문서키')
    b = b[~b['중복여부'] & ~b['용도군'].isin(a2.EXCLUDE_GROUPS)].sort_values('문서키').drop_duplicates(['제조사', '제품명'])
    b['순위가능'] = b['OELcov_TierB'] >= a2.MIN_COVERAGE
    b['석유류'] = b['Tfp_Liaw_C'].map(petroleum_class)
    r = b[b['순위가능']].dropna(subset=['Tfp_Liaw_C'])
    res = {'rankable_n': int(len(r)),
           'spearman_logI_vs_Tfp_all': round(float(np.log(r['I_mix_TierB']).corr(r['Tfp_Liaw_C'], method='spearman')), 3),
           'spearman_logI_vs_logF25_all': round(float(np.log(r['I_mix_TierB']).corr(np.log(r['F25']), method='spearman')), 3)}
    rows, quad_rows, dom_rows = [], [], []
    for g, d in r.groupby('비교군'):
        if len(d) < 3:
            continue
        rho = float(np.log(d['I_mix_TierB']).corr(d['Tfp_Liaw_C'], method='spearman'))
        rows.append({'비교군': g, 'n': len(d), 'spearman(logI, Tfp)': round(rho, 3)})
        mi, mt = d['I_mix_TierB'].median(), d['Tfp_Liaw_C'].median()
        for t in d.itertuples(index=False):
            low_h, low_f = t.I_mix_TierB <= mi, t.Tfp_Liaw_C >= mt
            quad_rows.append({'비교군': g, '제품명': t.제품명, '제조사': t.제조사, 'I_mix_TierB': round(t.I_mix_TierB, 1),
                              'Tfp_Liaw_C': round(t.Tfp_Liaw_C, 1), 'F25': round(t.F25, 3), '석유류': t.석유류,
                              '사분면': ('둘 다 낮음' if low_h and low_f else '노출만 낮음' if low_h else
                                      '화재만 낮음' if low_f else '둘 다 높음')})
        # 확실한 대체: 노출 P95 < 상대 P5, 인화점 P5 > 상대 P95
        for a_ in d.itertuples(index=False):
            cands = [b_ for b_ in d.itertuples(index=False)
                     if b_.I_mix_TierB_P95 < a_.I_mix_TierB_P5 and b_.Tfp_P5 > a_.Tfp_P95]
            if cands:
                best = max(cands, key=lambda c: (a_.I_mix_TierB / c.I_mix_TierB) * (1 + (c.Tfp_Liaw_C - a_.Tfp_Liaw_C) / 100))
                dom_rows.append({'비교군': g, '현재제품': a_.제품명, '현재 I': round(a_.I_mix_TierB, 1),
                                 '현재 Tfp': round(a_.Tfp_Liaw_C, 1), '현재 석유류': a_.석유류,
                                 '대체제품': best.제품명, '대체 I': round(best.I_mix_TierB, 1),
                                 '대체 Tfp': round(best.Tfp_Liaw_C, 1), '대체 석유류': best.석유류,
                                 'I 감소배수': round(a_.I_mix_TierB / best.I_mix_TierB, 2),
                                 'Tfp 상승_C': round(best.Tfp_Liaw_C - a_.Tfp_Liaw_C, 1),
                                 '석유류 구분 변경': a_.석유류 != best.석유류})
    pd.DataFrame(rows).to_csv(OUT / 'RQ7_within_group_correlation.csv', index=False, encoding='utf-8-sig')
    quad = pd.DataFrame(quad_rows)
    quad.to_csv(OUT / 'RQ7_quadrants.csv', index=False, encoding='utf-8-sig')
    dom = pd.DataFrame(dom_rows)
    dom.to_csv(OUT / 'RQ8_dominating_alternatives.csv', index=False, encoding='utf-8-sig')
    wg = pd.DataFrame(rows)
    res.update({'groups_with_>=3_rankable': int(len(wg)),
                'within_group_spearman_median': round(float(wg['spearman(logI, Tfp)'].median()), 3) if len(wg) else None,
                'quadrant_counts': quad['사분면'].value_counts().to_dict() if len(quad) else {},
                'products_with_dominating_alt': int(len(dom)),
                'dominating_I_reduction_median': round(float(dom['I 감소배수'].median()), 2) if len(dom) else None,
                'dominating_Tfp_gain_median_C': round(float(dom['Tfp 상승_C'].median()), 1) if len(dom) else None,
                'dominating_class_change_n': int(dom['석유류 구분 변경'].sum()) if len(dom) else 0,
                'petroleum_class_counts_all': b['석유류'].value_counts().to_dict(),
                'all_unique_products': int(len(b)),
                'unrankable_but_flammable': int(((~b['순위가능']) & (b['F25'] >= 1)).sum()),
                'unrankable_n': int((~b['순위가능']).sum())})
    # 노출지수는 낮아 보이는데 화재는 높은 경우: OEL 공백 제품 + 일반 제품 (사분면 '노출만 낮음' & 제1석유류)
    if len(quad):
        res['low_exposure_but_class1_n'] = int(((quad['사분면'] == '노출만 낮음') & (quad['석유류'] == '제1석유류')).sum())
    b.to_csv(OUT / 'RQ7_product_table.csv', index=False, encoding='utf-8-sig')
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s, prod, _sp = a2.load()
    summ = a2.summarise(s)
    f = pd.read_csv(RES / 'flammability_scenarios_270.csv.gz')
    fs = fire_summary(f)
    msds = pd.read_csv(V2 / 'data' / 'msds_flash.csv')
    tab, _ = rq6(fs, prod, msds)
    print(tab.to_string(index=False))
    res = rq7_rq8(fs, summ, prod)
    res.update(rq6_bounds(prod, msds, _sp))
    # 모델 간 인화점 차이
    piv = fs.pivot(index='문서키', columns='model', values='Tfp_Liaw_C')
    res['Tfp_HANNA_minus_Raoult_median_C'] = round(float((piv['HANNA'] - piv['Raoult']).median()), 2)
    res['Tfp_HANNA_minus_UNIFAC_median_C'] = round(float((piv['HANNA'] - piv['UNIFAC']).median()), 2)
    (OUT / 'headline_flammability.json').write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(res, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
