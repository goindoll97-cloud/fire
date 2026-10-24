# -*- coding: utf-8 -*-
"""본문 결과 1·2와 사업장용 제품 목록.

P1  용도를 구분하지 않고 모은 제품 전체에서 노출 지수와 인화점의 방향 (일반 결과)
P2  감사를 통과한 같은 비교군 안의 모든 제품쌍에서 두 지표가 같은 방향인 비율 (같은 용도 결과)
P3  비교군별 제품 목록: 다른 제품에 두 지표 모두 밀리는지(지배) 여부 (사업장용, 탐색적)

P1은 용도 비교가 아니므로 대체 후보를 뜻하지 않는다. P2·P3는 용도감사 `탐색비교허용=1`인 제품만 쓰며
TDS·건조속도·수지 호환성 검증 전이므로 대체 승인이 아니다.
"""
from __future__ import annotations

import itertools
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyze_core as ac
import analyze_flammability as af
import analyze_v2 as a2

V2, RES = a2.V2, a2.RES
OUT = RES / 'core'
N_BOOT = 2000
SEED = 20261024
NOT_THINNER = {'검색 오탐·비신나'}


def prepare(model='HANNA', flam_file='flammability_scenarios_270.csv.gz'):
    s, prod, sp = a2.load()
    summ = a2.summarise(s)
    fs = af.fire_summary(pd.read_csv(RES / flam_file))
    msds = pd.read_csv(V2 / 'data' / 'msds_flash.csv')[['기존제품ID', '인화점_C']].rename(columns={'인화점_C': 'FP_MSDS'})
    prod = ac.refined_groups(prod, sp).merge(msds, on='기존제품ID', how='left')
    b = ac.model_table(summ, fs, prod, model)
    b = b[(~b['중복여부']) & ~b['용도군'].isin(NOT_THINNER)].sort_values('문서키')
    b = b.drop_duplicates(['제조사', '제품명']).copy()
    b['순위가능'] = b['OELcov_TierB'] >= a2.MIN_COVERAGE
    b['석유류'] = b['Tfp_Liaw_C'].map(af.petroleum_class)
    return b


def spearman(x, y):
    return float(pd.Series(x).corr(pd.Series(y), method='spearman'))


def boot_ci(df, fn, rng):
    vals = []
    n = len(df)
    for _ in range(N_BOOT):
        v = fn(df.iloc[rng.integers(0, n, n)])
        if np.isfinite(v):
            vals.append(v)
    return [round(float(np.percentile(vals, 2.5)), 3), round(float(np.percentile(vals, 97.5)), 3)] if vals else [None, None]


def p1(bm, rng, fig=True):
    rows = []
    for label, d in (('전체', bm), *((f'제조사 {m}', g) for m, g in bm.groupby('제조사'))):
        for fire_col, name in (('Tfp_Liaw_C', '예측 인화점'), ('FP_MSDS', 'MSDS 인화점')):
            dd = d.dropna(subset=[fire_col])
            if len(dd) < 5:
                continue
            f = lambda x, c=fire_col: spearman(np.log(x['I_mix_TierB']), x[c])
            rows.append({'대상': label, '인화점 기준': name, 'n': len(dd), 'Spearman(logI, 인화점)': round(f(dd), 3),
                         '95% CI(부트스트랩)': boot_ci(dd, f, rng)})
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / 'P1_pooled_relation.csv', index=False, encoding='utf-8-sig')
    if fig:
        col = {'제1석유류': '#D55E00', '제2석유류': '#0072B2', '제3석유류 이상': '#009E73'}
        plt.rcParams.update({'font.family': ['Noto Sans CJK KR', 'Noto Sans CJK JP', 'DejaVu Sans'], 'axes.unicode_minus': False})
        fig_, ax = plt.subplots(figsize=(7.2, 4.8), dpi=200)
        for k, d in bm.dropna(subset=['Tfp_Liaw_C']).groupby('석유류'):
            ax.scatter(d['Tfp_Liaw_C'], d['I_mix_TierB'], s=22, alpha=.8, color=col.get(k, '#555'),
                       edgecolor='white', linewidth=.5, label=f'{k} ({len(d)})')
        ax.axvline(af.CLASS1, color='#444', lw=1, ls='--')
        ax.set_yscale('log')
        ax.set_xlabel('예측 인화점 (°C, HANNA + Liaw)')
        ax.set_ylabel('흡입 노출 지수 I_mix (노출기준 배수, 포화증기)')
        ax.set_title('용도를 구분하지 않아도 인화점이 높을수록 노출 지수가 낮다', fontsize=10, loc='left')
        for sp_ in ('top', 'right'):
            ax.spines[sp_].set_visible(False)
        ax.grid(alpha=.2)
        ax.legend(frameon=False, fontsize=8)
        fig_.tight_layout()
        fig_.savefig(OUT / 'fig_core2_pooled.png')
        plt.close(fig_)
    return tab


def pairs_table(bg):
    rows = []
    for g, d in bg.groupby('비교군'):
        for a, b in itertools.combinations(d.itertuples(index=False), 2):
            dI = np.log(a.I_mix_TierB / b.I_mix_TierB)
            dT = a.Tfp_Liaw_C - b.Tfp_Liaw_C
            if dI == 0 or dT == 0:
                continue
            dM = (a.FP_MSDS - b.FP_MSDS) if pd.notna(a.FP_MSDS) and pd.notna(b.FP_MSDS) else np.nan
            robust = ((a.I_mix_TierB_P95 < b.I_mix_TierB_P5 or b.I_mix_TierB_P95 < a.I_mix_TierB_P5) and
                      (a.Tfp_P95 < b.Tfp_P5 or b.Tfp_P95 < a.Tfp_P5))
            rows.append({'비교군': g, '제품A': a.제품명, '제품B': b.제품명, '같은 제조사': a.제조사 == b.제조사,
                         '일치(예측)': (dI > 0) != (dT > 0), '일치(MSDS)': np.nan if (np.isnan(dM) or dM == 0) else ((dI > 0) != (dM > 0)),
                         '두 지표 모두 확실히 분리': bool(robust)})
    return pd.DataFrame(rows)


def p2(bm, rng, save=True):
    bg = bm[bm['탐색비교허용'] & bm['순위가능']].dropna(subset=['Tfp_Liaw_C'])
    bg = bg[bg.groupby('비교군')['문서키'].transform('size') >= 2]
    pairs = pairs_table(bg)
    if save:
        pairs.to_csv(OUT / 'P2_pairs.csv', index=False, encoding='utf-8-sig')
    share = lambda d: float(pairs_table(d)['일치(예측)'].mean()) if len(d) > 1 else np.nan
    by_g = pairs.groupby('비교군').agg(쌍수=('일치(예측)', 'size'), 일치비율=('일치(예측)', 'mean')).reset_index()
    by_g['일치비율'] = by_g['일치비율'].round(3)
    if save:
        by_g.to_csv(OUT / 'P2_pairwise_by_group.csv', index=False, encoding='utf-8-sig')
    # 제품 단위 부트스트랩 (같은 군 안에서만 쌍을 만든다)
    vals = []
    ids = bg.index.to_numpy()
    for _ in range(N_BOOT // 4):
        samp = bg.loc[rng.choice(ids, len(ids), replace=True)].reset_index(drop=True)
        samp = samp.drop_duplicates() if False else samp
        v = pairs_table(samp)['일치(예측)'].mean() if len(samp) > 1 else np.nan
        if np.isfinite(v):
            vals.append(v)
    ci = [round(float(np.percentile(vals, 2.5)), 3), round(float(np.percentile(vals, 97.5)), 3)]
    pm = pairs['일치(MSDS)'].dropna().astype(bool)
    rob = pairs[pairs['두 지표 모두 확실히 분리']]
    res = {'products': int(len(bg)), 'groups': int(bg['비교군'].nunique()), 'pairs': int(len(pairs)),
           'concordance_model_fp': round(float(pairs['일치(예측)'].mean()), 3), 'concordance_model_fp_ci95': ci,
           'concordance_msds_fp': round(float(pm.mean()), 3), 'pairs_msds_evaluable': int(len(pm)),
           'pairs_robust_n': int(len(rob)), 'concordance_robust': round(float(rob['일치(예측)'].mean()), 3) if len(rob) else None,
           'concordance_same_maker': round(float(pairs[pairs['같은 제조사']]['일치(예측)'].mean()), 3),
           'concordance_diff_maker': round(float(pairs[~pairs['같은 제조사']]['일치(예측)'].mean()), 3),
           'pairs_same_maker_n': int(pairs['같은 제조사'].sum())}
    return res, bg


def p3(bg):
    rows = []
    for g, d in bg.groupby('비교군'):
        for a in d.itertuples(index=False):
            dominated = any(b.I_mix_TierB < a.I_mix_TierB and b.Tfp_Liaw_C > a.Tfp_Liaw_C for b in d.itertuples(index=False))
            dominates = any(b.I_mix_TierB > a.I_mix_TierB and b.Tfp_Liaw_C < a.Tfp_Liaw_C for b in d.itertuples(index=False))
            label = ('두 지표 모두 더 나은 제품이 있음' if dominated else
                     '이 군에서 두 지표 모두 우위' if dominates else '다른 제품과 일장일단 (상충)')
            rows.append({'비교군': g, '제품명': a.제품명, '제조사': a.제조사, '노출 지수 I': round(a.I_mix_TierB, 1),
                         '예측 인화점_C': round(a.Tfp_Liaw_C, 1), 'MSDS 인화점_C': a.FP_MSDS, '석유류(모형)': a.석유류,
                         '군 내 판정': label, '상태': '탐색적(대체 승인 아님)'})
    out = pd.DataFrame(rows).sort_values(['비교군', '노출 지수 I'])
    out.to_csv(OUT / 'P3_field_list.csv', index=False, encoding='utf-8-sig')
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    b = prepare('HANNA')
    bm = b[b['순위가능']].dropna(subset=['Tfp_Liaw_C'])
    t1 = p1(bm, rng)
    r2, bg = p2(b, rng)
    t3 = p3(bg)
    all_row = t1[(t1['대상'] == '전체')].set_index('인화점 기준')
    res = {'pooled_n': int(len(bm)), 'pooled_unique_products_total': int(len(b)),
           'pooled_spearman_model_fp': float(all_row.loc['예측 인화점', 'Spearman(logI, 인화점)']),
           'pooled_ci_model_fp': all_row.loc['예측 인화점', '95% CI(부트스트랩)'],
           'pooled_spearman_msds_fp': float(all_row.loc['MSDS 인화점', 'Spearman(logI, 인화점)']),
           'pooled_ci_msds_fp': all_row.loc['MSDS 인화점', '95% CI(부트스트랩)'],
           **{f'pairwise_{k}': v for k, v in r2.items()},
           'field_list_labels': t3['군 내 판정'].value_counts().to_dict()}
    (OUT / 'headline_pooled.json').write_text(json.dumps(res, ensure_ascii=False, indent=2, default=str), encoding='utf-8')
    print(t1.to_string(index=False))
    print(json.dumps(res, ensure_ascii=False, indent=2, default=str))


if __name__ == '__main__':
    main()
