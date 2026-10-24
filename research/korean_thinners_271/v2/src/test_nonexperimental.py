# -*- coding: utf-8 -*-
"""Model-free regression tests for reviewer-facing inferential claims."""
import unittest
import numpy as np
import pandas as pd

import analyze_pooled as ap
import audit_revision_sensitivity as ars


class TestNonexperimental(unittest.TestCase):
    def test_pareto_labels_are_mutually_exclusive(self):
        # A dominates B; C trades off against A and B.
        d = pd.DataFrame({
            '비교군': ['test'] * 3, '문서키': ['A', 'B', 'C'],
            '제품명': ['A', 'B', 'C'], '제조사': ['T'] * 3,
            'I_mix_TierB': [1., 3., 2.],
            'Tfp_Liaw_C': [10., 5., 12.],
            'FP_MSDS': [10., 5., 12.], '석유류': ['test'] * 3,
        })
        oldout = ap.OUT
        try:
            from tempfile import TemporaryDirectory
            from pathlib import Path
            with TemporaryDirectory() as temp:
                ap.OUT = Path(temp)
                result = ap.p3(d).set_index('제품명')
                self.assertFalse(bool(result.loc['B', 'Pareto_비지배']))
                self.assertTrue(bool(result.loc['A', 'Pareto_비지배']))
                self.assertTrue(bool(result.loc['C', 'Pareto_비지배']))
                self.assertEqual(int(result.loc['B', '동시개선_상대제품수']), 2)
                self.assertEqual(int(result.loc['A', '동시우위_상대제품수']), 1)
                self.assertEqual(int(result.loc['C', '동시우위_상대제품수']), 1)
        finally:
            ap.OUT = oldout

    def test_pair_count_not_independent_observations(self):
        d = pd.DataFrame({
            '비교군': ['one'] * 3, '문서키': list('abc'),
            '제품명': list('abc'), '제조사': ['m'] * 3,
            'I_mix_TierB': [1., 2., 3.],
            'I_mix_TierB_P5': [0.8, 1.8, 2.8],
            'I_mix_TierB_P95': [1.2, 2.2, 3.2],
            'Tfp_Liaw_C': [30., 20., 10.],
            'Tfp_P5': [29., 19., 9.],
            'Tfp_P95': [31., 21., 11.], 'FP_MSDS': [30., 20., 10.]
        })
        self.assertEqual(len(ap.pairs_table(d)), 3)
        self.assertTrue(ap.pairs_table(d)['일치(예측)'].all())

    def test_revision_name_normalization_only(self):
        self.assertEqual(ars.normalized_name('  CT - 101 일반 신나  '),
                         'ct-101일반신나')
        self.assertNotEqual(ars.normalized_name('아크릴 신나'),
                            ars.normalized_name('에폭시 신나'))


if __name__ == '__main__':
    unittest.main()
