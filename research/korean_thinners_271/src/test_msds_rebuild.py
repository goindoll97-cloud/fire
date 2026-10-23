# -*- coding: utf-8 -*-
import csv
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
import fitz
import build_msds_registry as base
import audit_research_readiness as check

class TestResearchDataIntegrity(unittest.TestCase):
    def test_cas_checksum(self):
        self.assertTrue(base.valid_cas('1330-20-7'))
        self.assertTrue(base.valid_cas('78-93-3'))
        self.assertFalse(base.valid_cas('1330-20-8'))

    def test_noru_korean_chemical_inventory_suffix(self):
        self.assertEqual(base.parse_interval('/ KE-12345\n5~15\nnext material')[:3], (5.0,15.0,'RANGE_EXTRACTED_UNVERIFIED'))
        self.assertEqual(base.parse_interval('/ KE-1200 0~8')[:2],(0.0,8.0))
        self.assertEqual(base.parse_interval('61∼71')[:2],(61.0,71.0))

    def test_samhwa_english_split_lines(self):
        self.assertEqual(base.parse_interval('More than\n0.1 ~ Under\n1 %')[:2],(.1,1.0))
        self.assertEqual(base.parse_interval('60 이상 ~ 70\n% 미만')[:2],(60.0,70.0))
        self.assertEqual(base.parse_interval('nonsense')[2],'REVIEW_CONCENTRATION')

    def test_mixed_isomer_blocked(self):
        self.assertIn('1330-20-7',base.MIXTURE_CAS)
        self.assertIn('64742-95-6',base.MIXTURE_CAS)
        self.assertNotIn('78-93-3',base.MIXTURE_CAS)

    def test_continued_section3_after_section4(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            doc=fitz.open()
            p=doc.new_page();p.insert_text((50,80), '1.CHEMICAL PRODUCT INFORMATION\nProduct Name: Test Thinner\n3.CHEMICAL COMPOSITION\nToluene\n108-88-3\nMore than 40 ~ Under 60 %\n4.HEALTH HAZARD INFORMATION',fontsize=10)
            p=doc.new_page();p.insert_text((50,80),'4.HEALTH HAZARD INFORMATION\n3.CHEMICAL COMPOSITION\nMEK\n78-93-3\nMore than 40 ~ Under 60 %\n5.FIRE FIGHTING MEASURES',fontsize=10)
            sam=tmp/'samhwa.zip';nor=tmp/'noru.zip'
            with zipfile.ZipFile(sam,'w') as z:z.writestr('sample.PDF',doc.tobytes())
            with zipfile.ZipFile(nor,'w') as z:pass
            docs,comps,evidence=base.index_pdfs({'NORU':nor,'SAMHWA':sam},tmp)
            self.assertEqual(len(docs),1)
            self.assertEqual({x['CAS'] for x in comps},{'108-88-3','78-93-3'})
            self.assertEqual(docs[0]['함량합계100가능_표기성분만'],'예')

    def test_property_approval_requires_source_and_no_mix(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/'prop.csv'
            rr=[{'CAS':'1330-20-7','혼합물·이성질체':'예','SMILES':'Cc1ccccc1','HANNA입력승인':'승인'}]
            base.write_csv(f,rr,list(rr[0]))
            with self.assertRaisesRegex(ValueError,'혼합물'):
                check.approved_property_records(f)
            rr=[{'CAS':'78-93-3','혼합물·이성질체':'아니오','SMILES':'CCC(C)=O','HANNA입력승인':'승인',
                 '구조검증상태':'미검증','물성검증상태':'미검증','분자량_g_mol':'72','25C_포화증기압_kPa':'10',
                 '구조공식출처':'dummy','물성출처와측정조건':'dummy','실측또는DB근거문헌':'dummy','온도_K':'298.15'}]
            base.write_csv(f,rr,list(rr[0]))
            with self.assertRaisesRegex(ValueError,'미승인'):
                check.approved_property_records(f)

if __name__=='__main__':unittest.main()