# -*- coding: utf-8 -*-
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
sys.path.insert(0,str(Path(__file__).parent))
import auto_enrich_270 as a


def mock_fetch(url,params):
    if '/name/' in url and url.endswith('/cids/JSON'):
        if '999-00-9' in url:return b'{"IdentifierList":{"CID":[21,22]}}'
        return b'{"IdentifierList":{"CID":[180]}}'
    if '/property/' in url:return b'{"PropertyTable":{"Properties":[{"CID":180,"SMILES":"CC=O","MolecularWeight":"44.05","InChIKey":"TEST"}]}}'
    if '/synonyms/' in url:return b'{"InformationList":{"Information":[{"Synonym":["75-07-0"]}]}}'
    if 'pug_view' in url:return json.dumps({'Record':{'Section':[{'TOCHeading':'Chemical and Physical Properties','Section':[{'TOCHeading':'Vapor Pressure','Information':[{'Value':{'StringWithMarkup':[{'String':'100 mmHg at 20 C'}]}}]}]}]}}).encode()
    if 'lawSearch' in url:
        if params.get('target')=='admbyl':return b'<Search><item>annex links review</item></Search>'
        return '<AdmRulSearch><admrul><행정규칙명>화학물질 및 물리적 인자의 노출기준</행정규칙명><행정규칙일련번호>1234</행정규칙일련번호></admrul></AdmRulSearch>'.encode()
    if 'lawService' in url:return '<AdmRulService><행정규칙명>화학물질 및 물리적 인자의 노출기준</행정규칙명><별표>75-07-0 acetaldehyde, 20 ppm, check source</별표></AdmRulService>'.encode()
    if 'getChemList' in url:return '<response><body><item><casNo>75-07-0</casNo></item></body></response>'.encode()
    raise AssertionError(url)

class AutoEnrichTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.work=Path(self.temp.name)
        self.review=self.work/'02_review';self.review.mkdir()
        self.props=[{'CAS':'75-07-0','복합물또는이성질체':'아니오','SMILES':'','분자량_g_mol':'','25C_포화증기압_kPa':'','온도_K':'298.15','구조공식출처':'','물성출처와측정조건':'','구조검증상태':'미검증','물성검증상태':'미검증','HANNA입력승인':'미승인'},
         {'CAS':'1330-20-7','복합물또는이성질체':'예','SMILES':'','분자량_g_mol':'','25C_포화증기압_kPa':'','온도_K':'298.15','구조공식출처':'','물성출처와측정조건':'','구조검증상태':'미검증','물성검증상태':'미검증','HANNA입력승인':'미승인'}]
        a.write(self.review/a.PROPS,self.props,list(self.props[0]))
        self.twa=[{'CAS':x['CAS'],'한국_TWA_ppm':'','기준시점':'','고시명번호':'','공식출처':'','공식원문대조':'미검증','검증상태':'미검증'} for x in self.props]
        a.write(self.review/a.TWA,self.twa,list(self.twa[0]))
        prod=[{'문서키':'D1','공정수지호환성검증':'미검증','비교군ID':'','근거문서':''}]
        a.write(self.review/a.PRODUCT,prod,list(prod[0]))

    def test_pubchem_fills_only_unapproved_structure_and_mw(self):
        summary=a.property_enrich(fetch=mock_fetch,work=self.work,sleep=0,limit=2)
        rows=a.read(self.review/a.PROPS)
        self.assertEqual(rows[0]['SMILES'],'CC=O')
        self.assertEqual(rows[0]['분자량_g_mol'],'44.05')
        self.assertEqual(rows[0]['25C_포화증기압_kPa'],'')
        self.assertEqual(rows[0]['HANNA입력승인'],'미승인')
        self.assertEqual(rows[1]['SMILES'],'')
        self.assertEqual(summary['mixed_blocked'],1)
        self.assertTrue((self.work/'02a_pubchem'/'PubChem_54CAS_구조분자량_자동후보.csv').exists())

    def test_already_approved_row_is_untouched(self):
        rows=a.read(self.review/a.PROPS)
        rows[0]['HANNA입력승인']='승인'
        a.write(self.review/a.PROPS,rows,list(rows[0]))
        a.property_enrich(fetch=mock_fetch,work=self.work,sleep=0,limit=1)
        self.assertEqual(a.read(self.review/a.PROPS)[0]['SMILES'],'')

    def test_pubchem_must_crosscheck_cas_synonyms(self):
        result=a.cid_lookup('999-00-9',fetch=mock_fetch)
        self.assertEqual(result['상태'],'여러_CID_동정필요')
        result=a.cid_lookup('100-41-4',fetch=mock_fetch)
        self.assertEqual(result['상태'],'CAS_동의어_미확인')

    def test_experimental_vp_annotation_never_becomes_25c_numeric(self):
        a.property_enrich(fetch=mock_fetch,work=self.work,sleep=0,limit=1)
        r=a.read(self.work/'02a_pubchem'/'PubChem_54CAS_구조분자량_자동후보.csv')[0]
        self.assertIn('100 mmHg at 20 C',r['VP_후보문구'])
        self.assertEqual(a.read(self.review/a.PROPS)[0]['25C_포화증기압_kPa'],'')

    def test_official_twa_never_approves_kosha_reference(self):
        with patch.object(a,'load_env'),patch.dict(os.environ,{'LAW_OC':'testsecret','KOSHA_MSDS_SERVICE_KEY':'othersecret'}):
            summary=a.official_twa_evidence(fetch=mock_fetch,work=self.work,sleep=0,limit=1)
        self.assertEqual(summary['approved_Korean_TWA_added'],0)
        self.assertEqual(a.read(self.review/a.TWA)[0]['한국_TWA_ppm'],'')
        self.assertEqual(a.read(self.review/a.TWA)[0]['검증상태'],'미검증')
        text=' '.join(x.read_text(encoding='utf-8') for x in (self.work/'02b_twa_evidence').glob('*.xml'))
        self.assertNotIn('testsecret',text)
        self.assertNotIn('othersecret',text)
        self.assertTrue((self.work/'02b_twa_evidence'/'법제처_현행노출기준_행정규칙원문.xml').exists())

    def test_only_direct_25c_pressure_candidate_is_convertible(self):
        parsed,state=a.direct_25c_vp_candidates(['10 mmHg at 25 C'])
        self.assertEqual(state,'VP_25C_직접조건후보_검증필요')
        self.assertAlmostEqual(parsed[0],1.3332236842,places=7)
        self.assertIsNone(a.direct_25c_vp_candidates(['10 mmHg at 20 C'])[0])
        self.assertIsNone(a.direct_25c_vp_candidates(['10 mmHg at 25 C','20 mmHg at 25 C'])[0])

    def test_scientific_notation_vapor_pressure_does_not_parse_exponent_only(self):
        parsed,state=a.direct_25c_vp_candidates(
          ['4.60X10+1 mm Hg at 25 °C /Extrapolated/'])
        self.assertEqual(state,'VP_25C_직접조건후보_검증필요')
        self.assertAlmostEqual(parsed[0],46*0.13332236842105263,places=6)
        self.assertIsNone(a.direct_25c_vp_candidates(['4.60X10+1 mm Hg at 20 °C'])[0])

    def test_previously_wrong_unapproved_vp_is_quarantined_and_reassessed(self):
        rows=a.read(self.review/a.PROPS)
        rows[0]['CAS']='142-82-5'
        rows[0]['25C_포화증기압_kPa']='0.13332237'
        rows[0]['물성출처와측정조건']=(
          'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/8900/JSON'
          ' ; original: 4.60X10+1 mm Hg at 25 °C /Extrapolated/'
          ' ; 자동후보/출처 원문 미승인')
        a.write(self.review/a.PROPS,rows,list(rows[0]))
        def heptane_fetch(url, params):
            if '/name/' in url and url.endswith('/cids/JSON'):
                return b'{"IdentifierList":{"CID":[8900]}}'
            if '/property/' in url:
                return b'{"PropertyTable":{"Properties":[{"CID":8900,"SMILES":"CCCCCCC","MolecularWeight":"100.2"}]}}'
            if '/synonyms/' in url:
                return b'{"InformationList":{"Information":[{"Synonym":["142-82-5"]}]}}'
            if 'pug_view' in url:
                return json.dumps({'Record':{'Section':[{
                  'TOCHeading':'Vapor Pressure','Information':[{'Value':{
                   'StringWithMarkup':[{'String':'4.60X10+1 mm Hg at 25 °C /Extrapolated/'}]}}]}]}}).encode()
            raise AssertionError(url)
        summary=a.property_enrich(fetch=heptane_fetch,work=self.work,limit=1,sleep=0)
        row=a.read(self.review/a.PROPS)[0]
        self.assertEqual(summary['VP_candidates_corrected'],1)
        self.assertAlmostEqual(float(row['25C_포화증기압_kPa']),46*0.13332236842105263,places=5)
        self.assertEqual(row['HANNA입력승인'],'미승인')
        audit=a.read(self.work/'02a_pubchem'/'자동VP후보_오류정정감사.csv')
        self.assertEqual(audit[0]['기존값_kPa'],'0.13332237')

    def test_non_single_identifiers_are_never_structurally_approved(self):
        rows=a.read(self.review/a.PROPS)
        rows[0]['CAS']='8052-41-3'
        rows[1]['CAS']='26588-80-7'
        rows[0]['복합물또는이성질체']='아니오'
        rows[1]['복합물또는이성질체']='아니오'
        a.write(self.review/a.PROPS,rows,list(rows[0]))
        summary=a.property_enrich(fetch=mock_fetch,work=self.work,limit=2,sleep=0)
        self.assertEqual(summary['mixed_blocked'],2)
        result=a.read(self.review/a.PROPS)
        self.assertTrue(all(r['복합물또는이성질체']=='예' for r in result))
        self.assertTrue(all(r['SMILES']=='' and r['HANNA입력승인']=='미승인' for r in result))

    def test_queue_does_not_create_approvals(self):
        summary=a.report_queue(work=self.work)
        self.assertEqual(summary['review_queue_count'],5)
        self.assertEqual(summary['HANNA_approved_CAS'],0)
        self.assertEqual(summary['Korean_TWA_approved_CAS'],0)

if __name__=='__main__':unittest.main()