# -*- coding: utf-8 -*-
import json
import tempfile
import unittest
from pathlib import Path
import pandas as pd
import apply_korean_msds_oel as s

# Representative table lines from the user-provided GitHub legal JSON checked
# by CAS. This is a format smoke-test, not an official-source integration test.
VALUES = {
'71-43-2': ('벤젠','Benzene',0.5),
'67-56-1':('메틸 알코올','Methanol',200),
'67-64-1':('아세톤','Acetone',500),
'123-86-4':('노말-초산 부틸','n-Butyl acetate',150),
'141-78-6':('초산 에틸','Ethyl acetate',400),
'108-88-3':('톨루엔','Toluene',50),
'108-10-1':('메틸 이소부틸 케톤','Methyl isobutyl ketone',50),
}

class KoreanScreenTests(unittest.TestCase):
    def test_twa_row_parser(self):
        row='│226 │벤젠 │Benzene │C6H6 │0.5 │- │2.5 │- │[71-43-2] 발암성 1A │'
        self.assertEqual(s.parse_twa_row([row], '71-43-2')['KR_TWA_ppm'], 0.5)
        with self.assertRaises(ValueError):
            s.parse_twa_row([row, row], '71-43-2')
    def test_existing_model_results_korean_reweight(self):
        path=s.DEFAULT_INPUT
        if not path.is_file():
            path = Path(__file__).resolve().parent / '03_component_predictions(1).csv'
        if not path.is_file():self.skipTest('model output not yet generated')
        d=pd.read_csv(path,dtype={'product_id':str,'target_CAS':str,'CAS':str})
        lines=[]
        for i,(cas,(kr,en,twa)) in enumerate(VALUES.items(),1):
            lines.append(f'│{i} │{kr} │{en} │FORM │{twa} │- │- │- │[{cas}] │')
        cw=s.build_crosswalk(d,{'별표':{'별표단위':[{'별표번호':'0001','별표내용':[lines]}]}})
        o=s.apply_oel(d,cw)
        case=o[(o.product_id=='1183746')&(o.target_CAS=='71-43-2') & (o.target_wt_pct==0.039)&(o.model_id=='hanna')&(o.CAS=='71-43-2')]
        self.assertEqual(len(case),1)
        self.assertAlmostEqual(float(case.iloc[0]['KR_OEL_weighted_share_pct']),14.6759679,places=5)
        self.assertEqual(int(case.iloc[0]['KR_screening_rank']),3)
        self.assertEqual(o.groupby(s.KEYS).ngroups,90)

    def test_annex_duplicate_number_selects_canonical_key(self):
        doc = {'별표': {'별표단위': [
            {'별표번호': '0001', '별표키': '000100',
             '별표제목': '화학물질의 노출기준', '별표내용': [['│1 │벤젠 │Benzene │C6H6 │0.5 │- │2.5 │- │[71-43-2] │']]},
            {'별표번호': '0001', '별표키': '000102',
             '별표제목': '삭제', '별표내용': '<별표 1-2> <삭제>'},
        ]}}
        a = s.find_annex(doc, '0001')
        self.assertEqual(a['별표키'], '000100')
        self.assertEqual(len(s.annex_lines(a)), 1)

    def test_cas_synonyms_identical_limits(self):
        a='│194 │메틸 이소부틸 케톤 │Methyl isobutyl ketone │CH3COCH2CH(CH3)2 │50 │- │75 │- │[108-10-1] 발암성 2 │'
        b='│714 │헥손 │Hexone │CH3COCH2CH(CH3)2 │50 │- │75 │- │[108-10-1] 발암성 2 │'
        row = s.parse_twa_row([a, b], '108-10-1')
        self.assertEqual(row['KR_TWA_ppm'], 50.)
        self.assertEqual(row['KR_OEL_listing_count'], 2)
        self.assertEqual(row['KR_name_en'], 'Methyl isobutyl ketone')
        with self.assertRaisesRegex(ValueError, 'Conflicting Korean'):
            s.parse_twa_row([a, b.replace('│50 │', '│100 │')], '108-10-1')

    def test_real_json_integration_if_available(self):
        # The local D:\fire download is the authoritative schema in actual use.
        # This test only skips when the caller does not have those JSON files.
        if not (s.OEL_JSON.is_file() and s.MSDS_JSON.is_file()):
            self.skipTest('local Korean law JSON not present in this environment')
        doc, m = s.regulatory_document(s.OEL_JSON)
        self.assertEqual(s.find_annex(doc, '0001')['별표키'], '000100')
        law, msds_meta = s.regulatory_document(s.MSDS_JSON)
        self.assertGreater(len(s.verify_msds_law(law)), 10)
        self.assertEqual(msds_meta['발령번호'], '2026-26')
        df = pd.read_csv(s.DEFAULT_INPUT, dtype={'CAS':str})
        xw = s.build_crosswalk(df, doc)
        self.assertEqual(len(xw), len(df['CAS'].unique()))
        row = xw.loc[xw.CAS == '108-10-1'].iloc[0]
        self.assertEqual(row.KR_OEL_listing_count, 2)
    def test_law_annex6_verification(self):
        fake={'별표':{'별표단위':[{'별표번호':'0006','별표내용':[['┃ 건강  │7. 발암성 │0.1%┃','┃ 건강 │1. 급성 독성 │1%┃']]}]},'조문내용':{'x':'제11조(작성원칙) ... 별표 6'}}
        self.assertEqual(len(s.verify_msds_law(fake)),2)
        fake['별표']['별표단위'][0]['별표내용'][0][0]='┃발암성│1%┃'
        with self.assertRaises(ValueError):s.verify_msds_law(fake)

if __name__=='__main__':unittest.main()
