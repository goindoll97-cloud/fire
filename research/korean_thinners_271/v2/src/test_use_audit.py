# -*- coding: utf-8 -*-
"""Deterministic checks for the 270-document MSDS use-classification audit.

These tests do not run HANNA or claim actual process compatibility.
Run: python -m unittest discover -s src -p "test_use_audit.py" -v
"""
import unittest

import pandas as pd

import analyze_v2 as analysis


class UseAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.audit = pd.read_csv(analysis.USE_AUDIT, dtype=str, encoding='utf-8-sig').fillna('')
        cls.raw = pd.read_csv(
            analysis.V2.parent / 'data' / 'source_rebuild_20261023' / 'MSDS_원본대조_270.csv',
            dtype=str, encoding='utf-8-sig'
        ).fillna('')
        cls.source = pd.read_csv(
            analysis.V2.parent / 'data' / '제품271건_용도재분류.csv',
            dtype=str, encoding='utf-8-sig'
        ).fillna('')

    def test_document_and_source_integrity(self):
        a, raw = self.audit, self.raw
        self.assertEqual(len(a), 270)
        self.assertEqual(len(raw), 270)
        self.assertTrue(a['문서키'].is_unique)
        self.assertTrue(raw['문서키'].is_unique)
        self.assertEqual(set(a['문서키']), set(raw['문서키']))
        joined = a[['문서키', '기존제품ID', '원본SHA256']].merge(
            raw[['문서키', '기존제품ID', '원본SHA256']], on='문서키',
            suffixes=('_audit', '_original'), validate='one_to_one'
        )
        self.assertTrue((joined['기존제품ID_audit'] == joined['기존제품ID_original']).all())
        self.assertTrue((joined['원본SHA256_audit'] == joined['원본SHA256_original']).all())
        self.assertTrue(set(raw['기존제품ID']).issubset(set(self.source['제품ID'])))

    def test_no_unapproved_field_substitution(self):
        self.assertTrue((self.audit['현장대체검증'] == '0').all())
        self.assertTrue(self.audit['탐색비교허용'].isin(['0', '1']).all())
        held = self.audit['용도감사상태'].str.startswith('검증유보')
        self.assertTrue((self.audit.loc[held, '탐색비교허용'] == '0').all())

    def test_conflict_items_held(self):
        byid = self.audit.set_index('기존제품ID')
        for pid in ['NR002', 'NR024', 'NR025', 'NR071', 'NR095', 'NR096', 'NR097', 'NR138']:
            self.assertEqual(byid.loc[pid, '탐색비교허용'], '0')
        for pid in ['NR027', 'NR036', 'NR115', 'NR133']:
            self.assertEqual(byid.loc[pid, '탐색비교허용'], '0')

    def test_section_one_reclassifications(self):
        byid = self.audit.set_index('기존제품ID')
        for pid in ['NR011', 'NR148', 'NR150', 'NR151', 'NR179']:
            self.assertEqual(byid.loc[pid, '기존용도군'], '명칭만으로 판별 불가')
            self.assertEqual(byid.loc[pid, '용도군_감사'], '산업·기재 지정용')
            self.assertTrue(byid.loc[pid, 'MSDS제1항권고용도'])
        self.assertIn('내화', byid.loc['NR116', '비교세부군'])

    def test_fine_grained_groups_do_not_mix_known_processes(self):
        byid = self.audit.set_index('기존제품ID')
        pairs = [('NR061', 'NR058'), ('NR171', 'NR169'), ('NR166', 'NR016'),
                 ('NR063', 'NR172'), ('NR155', 'NR174')]
        for p1, p2 in pairs:
            self.assertNotEqual(byid.loc[p1, '비교군_감사'], byid.loc[p2, '비교군_감사'])

    def test_pipeline_apply(self):
        source_by_id = self.source.set_index('제품ID')
        prod = self.raw[['문서키', '기존제품ID']].copy()
        prod['용도군'] = prod['기존제품ID'].map(source_by_id['1차 용도분류'])
        prod['세부분류'] = prod['기존제품ID'].map(source_by_id['세부분류'])
        prod['비교군'] = 'OLD_UNVERIFIED_GROUP'
        after = analysis.apply_use_audit(prod)
        self.assertEqual(len(after), 270)
        self.assertEqual(int(after['탐색비교허용'].sum()), 165)
        self.assertFalse(after['현장대체검증'].any())
        self.assertFalse((after['비교군'] == 'OLD_UNVERIFIED_GROUP').any())

    def test_fail_closed_on_missing_document(self):
        partial = self.raw[['문서키', '기존제품ID']].iloc[:-1].copy()
        with self.assertRaises(ValueError):
            analysis.apply_use_audit(partial.assign(용도군='미확인'))


if __name__ == '__main__':
    unittest.main()
