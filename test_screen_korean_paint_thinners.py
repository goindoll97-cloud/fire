import csv
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np

import screen_korean_paint_thinners as m


def sample_line(n, ko, en, twa, cas):
    return '│'+ '│'.join([str(n),ko,en,'C6H6',str(twa),'-','-','-',f'[{cas}]'])+'│'


class DomesticTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.base=Path(self.tmp.name)
        oel={'AdmRulService':{'행정규칙기본정보':{'발령번호':'2020-48','시행일자':'20200116'},
                '별표':{'별표단위':[{'별표키':'000100','별표제목':'화학물질의 노출기준','별표내용':[[
                    sample_line(1,'톨루엔','Toluene',50,'108-88-3'),
                    sample_line(2,'아세톤','Acetone',500,'67-64-1'),
                ]]}, {'별표키':'000102','별표제목':'삭제','별표내용':[]} ]}}}
        self.oel_path=self.base/'oel.json'
        self.oel_path.write_text(json.dumps(oel,ensure_ascii=False),encoding='utf-8')
        self.oel,_=m.parse_korean_twa(self.oel_path)
        self.props={'108-88-3':{'MW_g_mol':92.13,'VP_25C_kPa':3.79,'SMILES':'Cc1ccccc1','chemical':'toluene'},
                    '67-64-1':{'MW_g_mol':58.08,'VP_25C_kPa':30.60,'SMILES':'CC(C)=O','chemical':'acetone'}}
        self.a=self.row('A','108-88-3','60','60')
        self.b=self.row('A','67-64-1','40','40')

    def tearDown(self):
        self.tmp.cleanup()

    def row(self,pid,cas,lo,hi):
        return dict(product_id=pid,product_name=pid,manufacturer='Test Corp',market_country='KR',
                product_evidence_url='https://example.test/verified',msds_revision='2026-01-01',
                source_type='manufacturer_msds',use_group='Test common use',compatibility_group='Family 1',
                compatibility_verified='yes',compatibility_evidence_url='https://example.test/approved',
                CAS=cas,component_name=cas,wt_pct_min=lo,wt_pct_max=hi,msds_disclosure='yes')

    def test_domestic_parser_not_niosh(self):
        self.assertEqual(self.oel['108-88-3']['TWA_ppm'],50)
        self.assertEqual(self.oel['67-64-1']['TWA_ppm'],500)
        self.assertEqual(len(self.oel),2)

    def test_domestic_non_domestic_excluded(self):
        a=dict(self.a,market_country='US')
        self.assertEqual(m.check_product([a,self.b],self.props,self.oel)[0], 'PRODUCT_METADATA_INCONSISTENT')
        aa=[dict(x,market_country='US') for x in (self.a,self.b)]
        self.assertEqual(m.check_product(aa,self.props,self.oel)[0],'NOT_VERIFIED_KOREAN_MARKET')

    def test_disclosure_range_excluded(self):
        a=dict(self.a,wt_pct_min='55',wt_pct_max='60')
        self.assertEqual(m.check_product([a,self.b],self.props,self.oel)[0],'CONCENTRATION_RANGE_NOT_RANKED')

    def test_incomplete_total_excluded(self):
        b=dict(self.b,wt_pct_min='30',wt_pct_max='30')
        self.assertEqual(m.check_product([self.a,b],self.props,self.oel)[0],'INCOMPLETE_COMPOSITION_MASS_BALANCE')

    def test_missing_domestic_oel_excluded(self):
        oel=dict(self.oel);oel.pop('108-88-3')
        self.assertIn('MISSING_OR_AMBIGUOUS_KOREAN_TWA',m.check_product([self.a,self.b],self.props,oel)[0])

    def test_actual_equilibrium_not_workplace(self):
        smry,detail=m.screen_composition([('108-88-3',60),('67-64-1',40)],self.props,self.oel,np.ones(2))
        self.assertEqual(len(detail),2)
        self.assertAlmostEqual(sum(x['predicted_vapor_mol_pct'] for x in detail),100)
        self.assertAlmostEqual(detail[0]['equilibrium_screening_priority'],
            detail[0]['predicted_partial_pressure_kPa']/101.325*1e6/50)
        self.assertTrue(all(x['priority_interpretation'].startswith('SCREENING_ONLY') for x in detail))

    def test_invalid_gamma_rejected(self):
        with self.assertRaises(ValueError):
            m.screen_composition([('108-88-3',60),('67-64-1',40)],self.props,self.oel,[0,1])

    def test_run_no_automatic_substitution_without_hanna(self):
        rows=[self.a,self.b]
        for x in (self.a,self.b):
            rows.append(dict(x,product_id='B',product_name='B'))
        inputs=self.base/'inputs.csv'
        with inputs.open('w',newline='',encoding='utf-8') as f:
            wr=csv.DictWriter(f,fieldnames=m.COLUMNS)
            wr.writeheader();wr.writerows(rows)
        ref=self.base/'reference.csv'
        with ref.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=['CAS','chemical','MW_g_mol','VP_25C_kPa','SMILES','OEL_NIOSH_ppm'])
            writer.writeheader()
            for k,v in self.props.items():
                writer.writerow(dict(CAS=k,**v,OEL_NIOSH_ppm='99999'))
        fake=types.ModuleType('compare_thinner_raoult_unifac_hanna')
        class Stub:
            def __init__(self,*_): pass
            def predict(self,mid,cas,xs,temp,props):
                if mid=='hanna': raise RuntimeError('HANNA model absent')
                return np.ones(len(cas))
        fake.GammaModels=Stub
        with patch.dict(sys.modules,{'compare_thinner_raoult_unifac_hanna':fake}):
            result=m.run(inputs,self.base/'out',['hanna','raoult'],self.base,self.oel_path,ref)
        self.assertEqual(result['model_ok'],2)
        self.assertEqual(result['candidates'],0)
        self.assertEqual(result['total'],2)
        self.assertEqual(m.read_csv(self.base/'out'/'05_eligible_selection_candidates.csv'),[])
        self.assertFalse(json.loads((self.base/'out'/'06_run_metadata.json').read_text())['us_niosh_oel_used'])


if __name__ == '__main__':
    unittest.main()