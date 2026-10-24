# -*- coding: utf-8 -*-
"""Document-version selection audit: re-screen existing model outputs only.

Not a manufacturer-approved equivalence check. Same product name does not
establish that two formulations, coating systems or MSDS revisions are identical.
Unknown dates do not outrank confirmed dates. Original raw MSDS files are intact.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

import analyze_core as ac
import analyze_flammability as af
import analyze_pooled as ap
import analyze_v2 as a2

OUT = a2.RES / 'core'
BASE = a2.V2.parent


def normalized_name(value):
    return re.sub(r'\s+', '', str(value).casefold())


def source_versions(prod):
    src = pd.read_csv(BASE / 'data' / '제품271건_용도재분류.csv',
                      dtype=str, encoding='utf-8-sig').fillna('')
    cols = src[['제품ID', 'MSDS 개정일', 'MSDS 제1항 권고용도']].rename(
        columns={'제품ID': '기존제품ID'})
    out = prod.merge(cols, on='기존제품ID', how='left', validate='one_to_one')
    out['동명군'] = out['제조사'].astype(str) + '|' + out['제품명'].map(normalized_name)
    out['개정일_파싱'] = pd.to_datetime(out['MSDS 개정일'], errors='coerce')
    latest = out.sort_values(['동명군', '개정일_파싱', '기존제품ID'],
                             ascending=[True, False, True], na_position='last')
    selected = set(latest.drop_duplicates('동명군')['문서키'])
    out['최신개정_대표'] = out['문서키'].isin(selected)
    # Flag conflicting extracted uses within same product name, not just version date.
    uses = out.groupby('동명군')['MSDS 제1항 권고용도'].transform(
        lambda s: len(set(v.strip() for v in s if str(v).strip())) > 1)
    out['동명군_권고용도충돌'] = uses
    return out


def summary(label, frame):
    a = frame[frame['OELcov_TierB'].ge(0.9) & frame['I_mix_TierB'].gt(0)].dropna(
        subset=['Tfp_Liaw_C'])
    p1 = float(np.log(a['I_mix_TierB']).corr(a['Tfp_Liaw_C'], method='spearman'))
    same = a[a['탐색비교허용']]
    same = same[same.groupby('비교군')['문서키'].transform('size').ge(2)]
    pairs = ap.pairs_table(same) if len(same) else pd.DataFrame()
    return {
        'selection_rule': label, 'pooled_n': len(a),
        'pooled_spearman': round(p1, 3) if np.isfinite(p1) else None,
        'same_use_products': len(same), 'groups': same['비교군'].nunique(),
        'same_use_pairs': len(pairs),
        'same_direction_share': round(float(pairs['일치(예측)'].mean()), 3)
        if len(pairs) else None,
        'interpretation': 'same-name latest-version sensitivity, NOT approved formulation interchangeability',
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s, prod, sp = a2.load()
    summ = a2.summarise(s)
    fs = af.fire_summary(pd.read_csv(a2.RES / 'flammability_scenarios_270.csv.gz'))
    msds = pd.read_csv(a2.V2 / 'data' / 'msds_flash.csv')[
        ['기존제품ID', '인화점_C']].rename(columns={'인화점_C': 'FP_MSDS'})
    doc = ac.refined_groups(prod, sp).merge(
        msds, on='기존제품ID', how='left', validate='one_to_one')
    doc = source_versions(doc)
    doc[['문서키', '기존제품ID', '제조사', '제품명', 'MSDS 개정일',
         '동명군', '동명군_권고용도충돌', '조성지문', '중복여부',
         '최신개정_대표', '용도감사상태', '탐색비교허용']].to_csv(
             OUT / 'S3_document_version_audit.csv', index=False, encoding='utf-8-sig')
    all_products = ac.model_table(summ, fs, doc, 'HANNA')
    old = all_products[~all_products['중복여부']].sort_values('문서키')
    old = old.drop_duplicates(['제조사', '제품명'])
    newer = all_products[all_products['최신개정_대표']].copy()
    # A stricter option removes cross-name identical composition fingerprints
    # AFTER the most recent revision has been selected, for sensitivity only.
    latest_comp = newer.drop_duplicates('조성지문')
    rows = [summary('existing_fingerprint_first', old),
            summary('latest_same_name', newer),
            summary('latest_then_fingerprint', latest_comp)]
    pd.DataFrame(rows).to_csv(
        OUT / 'S3_revision_sensitivity.csv', index=False, encoding='utf-8-sig')
    print('NONEXPERIMENTAL_REVISION_SENSITIVITY ' +
          json.dumps({'versions': len(doc),
                      'same_name_extra_documents': int(len(doc) - doc['동명군'].nunique()),
                      'missing_revision_dates': int(doc['개정일_파싱'].isna().sum()),
                      'conflicting_use_rows': int(doc['동명군_권고용도충돌'].sum()),
                      'results': rows}, ensure_ascii=False))


if __name__ == '__main__':
    main()
