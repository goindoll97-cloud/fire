# -*- coding: utf-8 -*-
"""S5: model-input surrogate use and vapor OEL coverage threshold sensitivity.

Excluding UVCB-based surrogate compositions is NOT a validation of the remaining
predictions; it is a robustness/scope check on dependence on surrogate assumptions.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

import analyze_pooled as ap

OUT = ap.OUT


def screen(b, label, threshold=0.9, surrogate=None):
    sub = b[b['OELcov_TierB'].ge(threshold) & b['I_mix_TierB'].gt(0)].dropna(
        subset=['Tfp_Liaw_C']).copy()
    if surrogate is not None:
        flag = sub['UVCB대리물질사용'].astype(str).str.lower().eq('true')
        sub = sub[flag.eq(surrogate)].copy()
    rho = np.log(sub['I_mix_TierB']).corr(
        sub['Tfp_Liaw_C'], method='spearman') if len(sub) >= 3 else np.nan
    same = sub[sub['탐색비교허용']]
    same = same[same.groupby('비교군')['문서키'].transform('size') >= 2]
    pairs = ap.pairs_table(same) if len(same) else pd.DataFrame()
    return {
        'condition': label, 'OEL_coverage_min': threshold,
        'surrogate_filter': 'all' if surrogate is None else ('UVCB surrogate used' if surrogate else 'no UVCB surrogate'),
        'pooled_n': len(sub),
        'pooled_spearman': round(float(rho), 3) if np.isfinite(rho) else None,
        'same_use_products': len(same),
        'same_use_groups': same['비교군'].nunique(),
        'same_use_pairs': len(pairs),
        'same_direction_share': round(float(pairs['일치(예측)'].mean()), 3) if len(pairs) else None,
        'interpretation': 'restricted-set screening robustness only, not experimental validation'
    }


def main():
    b = ap.prepare('HANNA')
    cases = [
        screen(b, 'baseline_TierB_90pct', 0.9, None),
        screen(b, 'TierB_95pct', 0.95, None),
        screen(b, 'TierB_99pct', 0.99, None),
        screen(b, 'TierB_90pct_excluding_UVCB_surrogates', 0.9, False),
        screen(b, 'TierB_90pct_UVCB_surrogate_subset', 0.9, True),
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(cases).to_csv(
        OUT / 'S5_coverage_uvcb_sensitivity.csv', index=False, encoding='utf-8-sig')
    print('NONEXPERIMENTAL_S5_COVERAGE_UVCB ' + json.dumps(cases, ensure_ascii=False))


if __name__ == '__main__':
    main()
