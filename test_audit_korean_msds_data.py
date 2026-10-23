# -*- coding: utf-8 -*-
import unittest
from copy import deepcopy
from pathlib import Path
from audit_korean_msds_data import audit,load,ROOT,DEFAULT_PRODUCTS,DEFAULT_COMPONENTS

class TestKoreanMSDSAudit(unittest.TestCase):
    def setUp(self):
        self.p=load(DEFAULT_PRODUCTS)
        self.c=load(DEFAULT_COMPONENTS)
    def test_seven_products_39_component_rows(self):
        result=audit(self.p,self.c)
        self.assertEqual(len(result),7)
        self.assertEqual(sum(x[1] for x in result),39)
        self.assertTrue(all(x[4] for x in result))
    def test_distinct_noroo_gold_and_dr170l(self):
        ids={r['product_id']:r['msds_number'] for r in self.p}
        self.assertNotEqual(ids['KR_NOROO_DR170_GOLD'],ids['KR_NOROO_DR170L'])
    def test_component_removed_fails(self):
        with self.assertRaises(ValueError):
            audit(self.p,self.c[:-1])
    def test_modified_bound_fails(self):
        rows=deepcopy(self.c)
        rows[0]['wt_pct_min']='0'
        with self.assertRaises(ValueError):
            audit(self.p,rows)
    def test_mixed_xylene_label_omission_fails(self):
        rows=deepcopy(self.c)
        for c in rows:
            if c['CAS']=='1330-20-7':
                c['chemical_identity']='single'
        with self.assertRaises(ValueError):
            audit(self.p,rows)
    def test_original_pdf_hash_optional(self):
        root=ROOT/'local_msds_pdfs'
        if not root.is_dir():
            self.skipTest('원본 PDF 7개는 GitHub에 보관하지 않음; --pdf-dir 별도 확인')
        self.assertEqual(len(audit(self.p,self.c,root)),7)

if __name__=='__main__':
    unittest.main()
