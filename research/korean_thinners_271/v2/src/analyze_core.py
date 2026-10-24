# -*- coding: utf-8 -*-
"""본문(core) 분석: '사업주 이익(화재)과 근로자 건강(흡입)은 같은 방향인가?'

핵심 질문  같은 용도·같은 계절/건조 등급의 신나 안에서, 인화점이 높은(화재 위험이 낮은) 제품은
           흡입 노출 지수도 낮은가? 어떤 제품이 예외인가? 바꾸면 얼마나 좋아지는가?
AI 역할    HANNA 활동도계수 하나로 화재 지표(Liaw 인화점)와 노출 지표(I_mix)를 동시에 계산한다.
           C4 표는 기존 방법(Raoult, UNIFAC)이 계산하지 못하거나 결과가 달라지는 지점만 보여 준다.

출력: results/core/  (표 CSV, 그림 PNG, headline_core.json)
의존: run_v2.py, run_flammability.py 의 결과 파일, data/msds_flash.csv
"""
from __future__ import annotations

import itertools
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyze_flammability as af
import analyze_v2 as a2

V2, RES = a2.V2, a2.RES
OUT = RES / 'core'
RETARDER_FP = 60.0          # 성분 인화점 기준 (°C)
RETARDER_SHARE = 0.6        # 이 몰분율 이상이 고인화점 성분이면 지연제형
MODELS = ('HANNA', 'Raoult', 'UNIFAC')


def refined_groups(prod, sp):
    """비교군 = 용도군 + 계절/건조 변형. 건조속도 조절용은 조성으로 지연제형/혼합 지건형을 분리한다."""
    fl = pd.read_csv(V2 / 'data' / 'pure_flammability.csv').set_index('CAS')['FP_C']
    s = sp[(sp.model == 'Raoult') & (sp['T'] == '25C')].copy()
    s['high'] = s['CAS'].map(fl) >= RETARDER_FP
    share = s.groupby('문서키').apply(lambda d: d.loc[d['high'], 'x_liq'].sum(), include_groups=False)
    prod = prod.copy()
    prod['고인화점성분_몰분율'] = prod['문서키'].map(share)
    ret = (prod['용도군'] == '건조속도 조절용') & (prod['고인화점성분_몰분율'] >= RETARDER_SHARE)
    prod['비교군'] = np.where(ret, prod['비교군'] + ' (지연제형)',
                           np.where(prod['용도군'] == '건조속도 조절용', prod['비교군'] + ' (혼합 지건형)', prod['비교군']))
    return prod


def model_table(summ, fs, prod, model):
    h = summ[(summ['T'] == '25C') & (summ.model == model)][
        ['문서키', 'I_mix_TierB', 'I_mix_TierB_P5', 'I_mix_TierB_P95', 'OELcov_TierB']]
    f = fs[fs.model == model].drop(columns='model')
    return prod.merge(h, on='문서키').merge(f, on='문서키')


def universe(b):
    b = b[(~b['중복여부']) & b['탐색비교허용'] & ~b['용도군'].isin(a2.EXCLUDE_GROUPS)].sort_values('문서키')
    b = b.drop_duplicates(['제조사', '제품명']).copy()
    b['순위가능'] = b['OELcov_TierB'] >= a2.MIN_COVERAGE
    b['석유류'] = b['Tfp_Liaw_C'].map(af.petroleum_class)
    return b


def quadrant(d, fire_col):
    """Do not classify missing MSDS flash points using an imputed value."""
    valid = d[fire_col].notna() & d['I_mix_TierB'].notna()
    result = pd.Series(index=d.index, dtype='object')
    if not valid.any():
        return result
    subset = d.loc[valid]
    mi, mf = subset['I_mix_TierB'].median(), subset[fire_col].median()
    lo_h = subset['I_mix_TierB'] <= mi
    hi_f = subset[fire_col] >= mf
    result.loc[valid] = np.where(lo_h & hi_f, '둘 다 낮음',
                                 np.where(lo_h, '노출만 낮음',
                                          np.where(hi_f, '화재만 낮음', '둘 다 높음')))
    return result


def c1_c2(b):
    r = b[b['순위가능']].dropna(subset=['Tfp_Liaw_C']).copy()
    rows, quads = [], []
    for g, d in r.groupby('비교군'):
        if len(d) < 3:
            continue
        d = d.copy()
        rho = float(np.log(d['I_mix_TierB']).corr(d['Tfp_Liaw_C'], method='spearman'))
        dm = d.dropna(subset=['FP_MSDS'])
        rho_m = float(np.log(dm['I_mix_TierB']).corr(dm['FP_MSDS'], method='spearman')) if len(dm) >= 3 else np.nan
        rows.append({'비교군': g, 'n': len(d), 'Spearman(logI, 예측 인화점)': round(rho, 3),
                     'Spearman(logI, MSDS 인화점)': round(rho_m, 3)})
        d['사분면(예측 인화점)'] = quadrant(d, 'Tfp_Liaw_C')
        d['사분면(MSDS 인화점)'] = quadrant(d, 'FP_MSDS')
        quads.append(d)
    wg = pd.DataFrame(rows)
    if not quads:
        raise RuntimeError('감사된 동일 용도 비교군(3개 이상)이 없습니다. 용도감사 및 그룹을 검토하세요.')
    q = pd.concat(quads)
    wg.to_csv(OUT / 'C1_group_relation.csv', index=False, encoding='utf-8-sig')
    q[['비교군', '제조사', '제품명', 'I_mix_TierB', 'Tfp_Liaw_C', 'FP_MSDS', 'F25', '석유류',
       '사분면(예측 인화점)', '사분면(MSDS 인화점)']].round(3).to_csv(
        OUT / 'C2_quadrants.csv', index=False, encoding='utf-8-sig')
    return r, wg, q


def c3_swaps(q):
    rows = []
    for g, d in q.groupby('비교군'):
        for a in d.itertuples(index=False):
            cands = [c for c in d.itertuples(index=False)
                     if c.I_mix_TierB_P95 < a.I_mix_TierB_P5 and c.Tfp_P5 > a.Tfp_P95]
            if not cands:
                continue
            best = max(cands, key=lambda c: a.I_mix_TierB / c.I_mix_TierB)
            rows.append({'비교군': g, '현재제품': a.제품명, '현재 I': round(a.I_mix_TierB, 1),
                         '현재 인화점': round(a.Tfp_Liaw_C, 1), '현재 석유류': a.석유류,
                         '대체제품': best.제품명, '대체 I': round(best.I_mix_TierB, 1),
                         '대체 인화점': round(best.Tfp_Liaw_C, 1), '대체 석유류': best.석유류,
                         'I 감소배수': round(a.I_mix_TierB / best.I_mix_TierB, 2),
                         '인화점 상승_C': round(best.Tfp_Liaw_C - a.Tfp_Liaw_C, 1),
                         '석유류 구분 변경': a.석유류 != best.석유류})
    sw = pd.DataFrame(rows)
    sw.to_csv(OUT / 'C3_swap_alternatives.csv', index=False, encoding='utf-8-sig')
    return sw


def rank_flip_share(x, y, groups):
    flips = total = 0
    for _, idx in groups.groupby(groups).groups.items():
        idx = list(idx)
        for i, j in itertools.combinations(idx, 2):
            da, db = x[i] - x[j], y[i] - y[j]
            if da == 0 or db == 0:
                continue
            total += 1
            flips += (da > 0) != (db > 0)
    return flips, total


def c4_ai(summ, fs, prod, sp):
    tabs = {m: universe(model_table(summ, fs, prod, m)).set_index('문서키') for m in MODELS}
    # UNIFAC 계산 불가 제품은 summ/fs에 행이 없다 -> 모집단은 HANNA 기준
    base = tabs['HANNA']
    dmc = set(sp[(sp['CAS'] == '616-38-6')]['문서키'])
    rows = [{'지표': '계산 가능 제품 수', **{m: int(len(tabs[m])) for m in MODELS}},
            {'지표': 'DMC 함유 제품 중 계산 가능', **{m: int(sum(k in dmc for k in tabs[m].index)) for m in MODELS}}]
    for ref in ('Raoult', 'UNIFAC'):
        t = tabs[ref]
        common = base.index.intersection(t.index)
        cb, ct = base.loc[common], t.loc[common]
        cls = (cb['Tfp_Liaw_C'] < af.CLASS1) != (ct['Tfp_Liaw_C'] < af.CLASS1)
        rows.append({'지표': f'HANNA vs {ref}: 21°C 석유류 구분이 달라지는 제품', 'HANNA': int(cls.sum()),
                     ref: f'공통 {len(common)}종 중'})
        rk = common.intersection(base[base['순위가능']].index)
        g = base.loc[rk, '비교군']
        fh, th = rank_flip_share(base.loc[rk, 'I_mix_TierB'].to_dict(), t.loc[rk, 'I_mix_TierB'].to_dict(), g)
        ff, tf = rank_flip_share(base.loc[rk, 'Tfp_Liaw_C'].fillna(0).to_dict(), t.loc[rk, 'Tfp_Liaw_C'].fillna(0).to_dict(), g)
        rows.append({'지표': f'HANNA vs {ref}: 같은 비교군 내 노출 순위가 뒤바뀐 쌍', 'HANNA': fh, ref: f'{th}쌍 중'})
        rows.append({'지표': f'HANNA vs {ref}: 같은 비교군 내 인화점 순위가 뒤바뀐 쌍', 'HANNA': ff, ref: f'{tf}쌍 중'})
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / 'C4_ai_necessity.csv', index=False, encoding='utf-8-sig')
    return tab


def figures(r, q):
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': ['Noto Sans CJK KR', 'Noto Sans CJK JP', 'DejaVu Sans'], 'axes.unicode_minus': False})
    colors = {'둘 다 낮음': '#0072B2', '둘 다 높음': '#6b6b6b', '노출만 낮음': '#E69F00', '화재만 낮음': '#CC79A7'}
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=200)
    for name, d in q.groupby('사분면(예측 인화점)'):
        ax.scatter(d['Tfp_Liaw_C'], d['I_mix_TierB'], s=26, alpha=.85, color=colors[name], edgecolor='white',
                   linewidth=.6, label=f'{name} ({len(d)})')
    ax.axvline(af.CLASS1, color='#444', lw=1, ls='--')
    ax.text(af.CLASS1 + 1, ax.get_ylim()[1] * .55, '21°C\n(제1석유류↔제2석유류)', fontsize=8, color='#444')
    ax.set_yscale('log')
    ax.set_xlabel('예측 인화점 (°C, HANNA + Liaw)')
    ax.set_ylabel('흡입 노출 지수 I_mix (포화증기 기준, 노출기준 배수)')
    ax.set_title('같은 용도 안에서 화재 위험과 흡입 노출은 대체로 같은 방향이다', fontsize=10, loc='left')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=.2)
    ax.legend(frameon=False, fontsize=8, title='비교군 중앙값 기준', title_fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'fig_core1_health_vs_fire.png')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s, prod, sp = a2.load()
    summ = a2.summarise(s)
    fs = af.fire_summary(pd.read_csv(RES / 'flammability_scenarios_270.csv.gz'))
    msds = pd.read_csv(V2 / 'data' / 'msds_flash.csv')[['기존제품ID', '인화점_C']].rename(columns={'인화점_C': 'FP_MSDS'})
    prod = refined_groups(prod, sp).merge(msds, on='기존제품ID', how='left')
    b = universe(model_table(summ, fs, prod, 'HANNA'))
    r, wg, q = c1_c2(b)
    pd.DataFrame({
        '용도군': prod['용도군'], '비교세부군': prod['비교세부군'],
        '용도감사상태': prod['용도감사상태'], '탐색비교허용': prod['탐색비교허용']
    }).groupby(['용도군', '비교세부군', '용도감사상태', '탐색비교허용'], dropna=False).size().rename('문서수').reset_index().to_csv(
        OUT / 'use_audit_comparison_group_counts.csv', index=False, encoding='utf-8-sig')
    sw = c3_swaps(q)
    ai = c4_ai(summ, fs, prod, sp)
    figures(r, q)
    b.to_csv(OUT / 'core_product_table_HANNA.csv', index=False, encoding='utf-8-sig')
    same = (q['사분면(예측 인화점)'].isin(['둘 다 낮음', '둘 다 높음'])).mean()
    msds_quadrants = q['사분면(MSDS 인화점)'].dropna()
    same_m = msds_quadrants.isin(['둘 다 낮음', '둘 다 높음']).mean() if len(msds_quadrants) else float('nan')
    res = {'products_in_universe': int(len(b)), 'rankable_in_groups_ge3': int(len(q)),
           'groups_ge3': int(len(wg)),
           'spearman_median_model_fp': round(float(wg['Spearman(logI, 예측 인화점)'].median()), 3),
           'spearman_median_msds_fp': round(float(wg['Spearman(logI, MSDS 인화점)'].median()), 3),
           'share_same_direction_model_fp': round(float(same), 3),
           'share_same_direction_msds_fp': round(float(same_m), 3),
           'n_msds_quadrant_evaluable': int(len(msds_quadrants)),
           'n_model_quadrant_evaluable': int(q['사분면(예측 인화점)'].notna().sum()),
           'use_audit_excluded_documents': int((~prod['탐색비교허용']).sum()),
           'use_audit_policy': 'exploratory_only_not_confirmed_substitution',
           'quadrants_model_fp': q['사분면(예측 인화점)'].value_counts().to_dict(),
           'quadrants_msds_fp': q['사분면(MSDS 인화점)'].value_counts().to_dict(),
           'swap_products': int(len(sw)),
           'swap_I_reduction_median': round(float(sw['I 감소배수'].median()), 2) if len(sw) else None,
           'swap_fp_gain_median_C': round(float(sw['인화점 상승_C'].median()), 1) if len(sw) else None,
           'swap_class_change': int(sw['석유류 구분 변경'].sum()) if len(sw) else 0,
           'unrankable_n': int((~b['순위가능']).sum()),
           'unrankable_F25_ge_1': int(((~b['순위가능']) & (b['F25'] >= 1)).sum()),
           'class_counts': b['석유류'].value_counts().to_dict()}
    (OUT / 'headline_core.json').write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(res, ensure_ascii=False, indent=2))
    print(ai.to_string(index=False))
    print(wg.to_string(index=False))


if __name__ == '__main__':
    main()
