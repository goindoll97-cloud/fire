# -*- coding: utf-8 -*-
"""Inventory ALL source CAS and explicitly show incompleteness of legal screening.

A missing rule in the verified-regulation register is never treated as evidence
that the substance or product is not regulated. Official NICS 2026-5 Annex full
XLSX still requires ingest and independent row-level effective-date review.
"""
from __future__ import annotations
import json
import pandas as pd
import analyze_v2 as a2

V2=a2.V2
OUT=a2.RES/'regulatory_gap'
COMP=V2.parent/'data'/'source_rebuild_20261023'/'MSDS_구성성분_1454.csv'
RULES=V2/'data'/'verified_regulatory_thresholds_20261024.csv'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    src=pd.read_csv(COMP,encoding='utf-8-sig',dtype={'CAS':str})
    doc=pd.read_csv(V2/'results'/'products_270.csv',encoding='utf-8-sig')
    rules=pd.read_csv(RULES,encoding='utf-8-sig',dtype={'CAS':str})
    source=src.groupby('CAS').agg(source_MSDS_rows=('문서키','size'),
                                  source_documents=('문서키','nunique'),
                                  source_lower_min=('함량하한_pct','min'),
                                  source_upper_max=('함량상한_pct','max')
                                  ).reset_index()
    rule_counts=rules.groupby('CAS').agg(
        verified_clause_count=('regulation_category','size'),
        verified_categories=('regulation_category',
                             lambda s: '|'.join(sorted(set(s)))),
        source_links=('source_url',
                      lambda s:'|'.join(sorted(set(s))))
    ).reset_index()
    audit=source.merge(rule_counts,on='CAS',how='left',validate='one_to_one')
    audit['cas_regulation_verification_status']=audit['verified_clause_count'].apply(
        lambda x: 'selected_CAS_clause_verified' if pd.notna(x)
        else 'NOT_YET_REVIEWED_DO_NOT_INFER_NONREGULATION')
    audit['legal_non_applicability_approved']=False
    audit['needs_full_current_annex_review']=True
    audit= audit.sort_values('source_documents',ascending=False)
    audit.to_csv(OUT/'L1_all_source_CAS_legal_coverage_ledger.csv',
                 index=False,encoding='utf-8-sig')
    known=set(rules.CAS)
    perdoc=src.groupby('문서키').agg(
        source_CAS_count=('CAS','nunique'),
        known_legal_reviewed_CAS_count=('CAS',lambda s:len(set(s)&known)),
        not_reviewed_CAS_count=('CAS',lambda s:len(set(s)-known)),
        not_reviewed_CAS_list=('CAS',lambda s:'|'.join(sorted(set(s)-known)))
    ).reset_index().merge(doc[['문서키','기존제품ID','제품명','제조사']],
                          on='문서키',validate='one_to_one')
    perdoc['complete_legal_applicability_established']=False
    perdoc.to_csv(OUT/'L2_270_MSDS_unreviewed_CAS.csv',
                  index=False,encoding='utf-8-sig')
    summary={
        'source_documents':len(perdoc),
        'source_CAS_total':len(source),
        'source_CAS_with_at_least_one_verified_clause':int(
            source.CAS.isin(known).sum()),
        'source_CAS_not_yet_reviewed':int((~source.CAS.isin(known)).sum()),
        'source_MSDS_with_any_nonreviewed_CAS':int(
            perdoc.not_reviewed_CAS_count.gt(0).sum()),
        'complete_regulation_clearances':0,
        'verified_rule_rows_including_CAS_not_in_source':len(rules),
        'full_annex_2026_5_status':'Official XLSX found, but full row-level annex not ingested or verified',
        'future_effective_2026_16_status':'Not in force on 2026-10-24; do not merge into present baseline'
    }
    (OUT/'L0_legal_scope_integrity.json').write_text(
        json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print('ALL_SOURCE_CAS_LEGAL_AUDIT '+json.dumps(summary,ensure_ascii=False))

if __name__=='__main__':
    main()
