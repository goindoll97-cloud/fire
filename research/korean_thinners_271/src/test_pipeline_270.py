# -*- coding: utf-8 -*-
"""Synthetic unit tests, no network, HANNA checkpoint or true approved CAS data."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).parent))
import pipeline_270_core as p

class CoreTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.work=self.root/'out'
        self.source=self.root/'repo_data'
        self.source.mkdir()
        self.review=self.work/'02_review'
        self.review.mkdir(parents=True)
        self.patch=patch.multiple(p,ROOT=self.root,WORK=self.work,SOURCE=self.source)
        self.patch.start();self.addCleanup(self.patch.stop)
        self.comp=[{'문서키':'A','기존제품ID':'X','CAS':'A1', 'CAS분류':'순물질후보_미검증',
                    '함량하한_pct':'30','함량상한_pct':'50','함량추출상태':'RANGE_EXTRACTED_UNVERIFIED'},
                   {'문서키':'A','기존제품ID':'X','CAS':'A2','CAS분류':'순물질후보_미검증',
                    '함량하한_pct':'50','함량상한_pct':'70','함량추출상태':'RANGE_EXTRACTED_UNVERIFIED'}]
        self.reviews=[{'문서키':'A','기존제품ID':'X','제조사':'S','제품명':'검증된 가상A','용도ID':'도료계',
                       '비교군ID':'compat','MSDS제1항용도검증':'승인',
                       'MSDS제3항함량검증':'승인','전체모델성분확정':'승인','중복개정판검토':'승인',
                       '독립제품검증':'승인','공정수지호환성검증':'승인',
                       'TDS성능제약검증':'승인','근거문서':'fixture'}]
        p.write(self.review/'제품_원문TDS_검증.csv',self.reviews,p.PRODUCT_FIELDS)
        self.props=[{'CAS':cas,'복합물또는이성질체':'아니오','SMILES':'CC', '분자량_g_mol':'100',
          '25C_포화증기압_kPa':'4','온도_K':'298.15','구조공식출처':'fixture','물성출처와측정조건':'fixture',
          '구조검증상태':'승인','물성검증상태':'승인','HANNA입력승인':'승인'} for cas in ['A1','A2']]
        p.write(self.review/'CAS_물성_검증.csv',self.props,p.PROPERTY_FIELDS)
        self.twa=[{'CAS':cas,'한국_TWA_ppm':'50','기준시점':'2026-10-23', '고시명번호':'fixture',
                 '공식출처':'fixture','공식원문대조':'승인','검증상태':'승인'} for cas in ['A1','A2']]
        p.write(self.review/'CAS_한국TWA_검증.csv',self.twa,p.TWA_FIELDS)

    def test_env_load_does_not_overwrite_existing(self):
        (self.root/'.env.example').write_text('LAW_OC=\nKOSHA_MSDS_SERVICE_KEY=\n')
        (self.root/'.env').write_text('LAW_OC=FAKEKEY\n')
        with patch.dict('os.environ',{'LAW_OC':'shellValue'}):
            p.stage00(None)
            self.assertEqual(p.os.environ['LAW_OC'],'shellValue')
        self.assertIn('FAKEKEY',(self.root/'.env').read_text())
        self.assertNotIn('FAKEKEY',(self.work/'00_environment'/'status.json').read_text())

    def test_no_automatic_approval(self):
        self.props[0]['HANNA입력승인']='미승인'
        p.write(self.review/'CAS_물성_검증.csv',self.props,p.PROPERTY_FIELDS)
        self.assertEqual(set(p.approved_properties()),{'A2'})

    def test_mixture_rejected(self):
        self.props[0]['복합물또는이성질체']='예'
        p.write(self.review/'CAS_물성_검증.csv',self.props,p.PROPERTY_FIELDS)
        with self.assertRaisesRegex(ValueError,'혼합물'):p.approved_properties()

    def test_unverified_twa_rejected(self):
        self.twa[0]['공식원문대조']='미검증'
        p.write(self.review/'CAS_한국TWA_검증.csv',self.twa,p.TWA_FIELDS)
        with self.assertRaisesRegex(ValueError,'원문증거'):p.approved_twa()

    def test_twa_missing_blocks_entire_scenario(self):
        self.fake_source()
        self.twa[1]['검증상태']='미검증'
        p.write(self.review/'CAS_한국TWA_검증.csv',self.twa,p.TWA_FIELDS)
        o=self.work/'05_predictions';o.mkdir(parents=True)
        p.write(o/'HANNA_Raoult_시나리오.csv',
           [{'문서키':'A','시나리오ID':1,'Raoult_총압_kPa':1.0,'HANNA_총압_kPa':2.0}],
           ['문서키','시나리오ID','Raoult_총압_kPa','HANNA_총압_kPa'])
        p.write(o/'HANNA_Raoult_성분.csv',
           [{'문서키':'A','시나리오ID':1,'CAS':'A1','HANNA_분압_kPa':1,'Raoult_분압_kPa':0.5},
            {'문서키':'A','시나리오ID':1,'CAS':'A2','HANNA_분압_kPa':1,'Raoult_분압_kPa':0.5}],
           ['문서키','시나리오ID','CAS','HANNA_분압_kPa','Raoult_분압_kPa'])
        p.jsonfile(o/'prediction_provenance.json',{'source_sha256':{x:p.sha(self.source/x) for x in(p.DOCFILE,p.COMPFILE)},
          'property_review_sha256':p.sha(self.review/'CAS_물성_검증.csv'),
          'product_review_sha256':p.sha(self.review/'제품_원문TDS_검증.csv')})
        p.stage06(None)
        self.assertEqual(p.read(self.work/'06_hazard'/'성분별_잠재지표_시나리오.csv'),[])
        self.assertEqual(len(p.read(self.work/'06_hazard'/'보류.csv')),1)

    def test_same_use_requires_two_independent_verified_products(self):
        # Single approved doc may NOT get ranked.
        self.fake_source()
        path=self.work/'06_hazard'/'성분별_잠재지표_시나리오.csv'
        p.write(path,[{'문서키':'A','시나리오ID':'1','HANNA_최대성분정규화평형분압':'10'}],
                ['문서키','시나리오ID','HANNA_최대성분정규화평형분압'])
        p.jsonfile(self.work/'06_hazard'/'hazard_provenance.json',
                   {'verified_twa_review_sha256':p.sha(self.review/'CAS_한국TWA_검증.csv')})
        # Synthetic fixtures bypass the production 270-document invariant only here.
        with patch.object(p,'ensure_source',return_value=(
             [{'문서키':'A','기존제품ID':'X','제품명_파일명기반':'가상A',
               '함량합계100가능_표기성분만':'예'}],self.comp,[])):
            p.jsonfile(self.work/'05_predictions'/'prediction_provenance.json',
                 {'product_review_sha256':p.sha(self.review/'제품_원문TDS_검증.csv')})
            p.stage07(None)
        self.assertEqual(p.read(self.work/'07_same_use'/'같은용도_낮은잠재지표_후보.csv'),[])

    def fake_source(self):
        p.write(self.source/p.DOCFILE,[{'문서키':'A','기존제품ID':'X','원본SHA256':'x',
        '제품명_파일명기반':'검증된 가상A','함량합계100가능_표기성분만':'예'}],
        ['문서키','기존제품ID','원본SHA256','제품명_파일명기반','함량합계100가능_표기성분만'])
        p.write(self.source/p.COMPFILE,self.comp,list(self.comp[0]))
        p.write(self.source/p.MASTERFILE,[{'CAS':'A1'},{'CAS':'A2'}],['CAS'])

    def test_source_verification_strict(self):
        self.fake_source()
        with self.assertRaisesRegex(ValueError,'270/1454/54'):
            p.ensure_source()

if __name__=='__main__':unittest.main()