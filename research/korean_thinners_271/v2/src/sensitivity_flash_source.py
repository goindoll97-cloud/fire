# -*- coding: utf-8 -*-
"""성분 인화점 출처 등급(B·C)의 불확실성이 결과를 바꾸는지 확인하는 민감도 분석 (Raoult 기준, HANNA 없이 재현 가능).

선행 실행:  python run_flammability.py --fp-offset {0,15,-15} --out-suffix {_sens0,_sensp15,_sensm15}
출력: results/core/S1_flash_source_sensitivity.csv
"""
import numpy as np
import pandas as pd

import analyze_flammability as af
import analyze_pooled as ap

V = {'기준(오프셋 0)': '_sens0', '등급 B·C 성분 +15°C': '_sensp15', '등급 B·C 성분 −15°C': '_sensm15'}


def main():
    rng = np.random.default_rng(ap.SEED)
    rows, base = [], None
    for label, suf in V.items():
        b = ap.prepare('Raoult', f'flammability_scenarios_270{suf}.csv.gz')
        bm = b[b['순위가능']].dropna(subset=['Tfp_Liaw_C'])
        rho = float(np.log(bm['I_mix_TierB']).corr(bm['Tfp_Liaw_C'], method='spearman'))
        r2, _ = ap.p2(b, rng, save=False)
        cls = b.set_index('문서키')['석유류']
        tf = b.set_index('문서키')['Tfp_Liaw_C']
        if base is None:
            base = (cls, tf)
        common = cls.index.intersection(base[0].index)
        rows.append({'조건': label, 'P1 Spearman(전체)': round(rho, 3), 'n': len(bm),
                     'P2 제품쌍 일치율': r2['concordance_model_fp'], 'P2 쌍수': r2['pairs'],
                     '인화점 중앙 변화_C(기준 대비)': round(float((tf[common] - base[1][common]).abs().median()), 2),
                     '인화점 최대 변화_C': round(float((tf[common] - base[1][common]).abs().max()), 1),
                     '석유류 구분이 바뀐 제품(기준 대비)': int((cls[common] != base[0][common]).sum())})
    out = pd.DataFrame(rows)
    out.to_csv(ap.OUT / 'S1_flash_source_sensitivity.csv', index=False, encoding='utf-8-sig')
    print(out.to_string(index=False))


if __name__ == '__main__':
    main()
