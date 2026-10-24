# -*- coding: utf-8 -*-
"""No-laboratory sensitivity analysis using existing HANNA/Raoult/UNIFAC runs.

All outputs are *screening* statistics, not measured human exposure or flash points.
This module never replaces missing OELs with an estimated occupational limit.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

import analyze_v2 as a2
import analyze_pooled as ap

OUT = a2.RES / 'core'


def model_tier_table():
    s, prod, _ = a2.load()
    summ = a2.summarise(s)
    results = []
    for model in ('HANNA', 'Raoult', 'UNIFAC'):
        base = ap.prepare(model)
        for tier in ('A', 'B'):
            rows = base.copy()
            if tier == 'A':
                health = summ[(summ['model'] == model) & (summ['T'] == '25C')][[
                    '문서키', 'I_mix_TierA', 'I_mix_TierA_P5', 'I_mix_TierA_P95', 'OELcov_TierA']]
                rows = rows.drop(columns=[
                    'I_mix_TierB', 'I_mix_TierB_P5', 'I_mix_TierB_P95', 'OELcov_TierB', '순위가능'])
                health = health.rename(columns={
                    'I_mix_TierA': 'I_mix_TierB', 'I_mix_TierA_P5': 'I_mix_TierB_P5',
                    'I_mix_TierA_P95': 'I_mix_TierB_P95', 'OELcov_TierA': 'OELcov_TierB'})
                rows = rows.merge(health, on='문서키', how='left', validate='one_to_one')
            rows['순위가능'] = (rows['OELcov_TierB'] >= 0.9) & (rows['I_mix_TierB'] > 0)
            filtered = rows[rows['순위가능']].dropna(subset=['Tfp_Liaw_C'])
            rho = float(np.log(filtered['I_mix_TierB']).corr(
                filtered['Tfp_Liaw_C'], method='spearman')) if len(filtered) >= 3 else np.nan
            group_rows = rows[rows['탐색비교허용'] & rows['순위가능']].dropna(subset=['Tfp_Liaw_C'])
            group_rows = group_rows[group_rows.groupby('비교군')['문서키'].transform('size') >= 2]
            pairtab = ap.pairs_table(group_rows) if len(group_rows) else pd.DataFrame()
            results.append({
                'model': model, 'OEL_Tier': tier,
                'pooled_n': len(filtered), 'pooled_spearman': round(rho, 3) if np.isfinite(rho) else None,
                'same_use_n': len(group_rows), 'same_use_groups': group_rows['비교군'].nunique(),
                'same_use_pairs': len(pairtab),
                'same_direction_share': round(float(pairtab['일치(예측)'].mean()), 3)
                if len(pairtab) else None,
                'coverage_rule': '>=90% saturated vapor with assigned OEL',
                'scope': 'model-and-OEL sensitivity, no real exposure measurement',
            })
    return pd.DataFrame(results)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = model_tier_table()
    results.to_csv(OUT / 'S2_models_oel_sensitivity.csv', index=False, encoding='utf-8-sig')
    print('NONEXPERIMENTAL_OEL_MODEL_SENSITIVITY ' +
          json.dumps(results.to_dict(orient='records'), ensure_ascii=False))


if __name__ == '__main__':
    main()
