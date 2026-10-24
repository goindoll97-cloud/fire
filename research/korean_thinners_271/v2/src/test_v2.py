# -*- coding: utf-8 -*-
"""Unit tests for the v2 pipeline (run: python -m unittest discover -s v2/src -p 'test_*.py')."""
import csv
import random
import sys
import unittest
from collections import defaultdict
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parent
sys.path.insert(0, str(SRC))
V2 = SRC.parent

import parse_korean_oel  # noqa: E402
import run_v2  # noqa: E402


def rows(p):
    with p.open(encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


class TestOELParse(unittest.TestCase):
    def test_known_values(self):
        idx = {}
        for r in parse_korean_oel.parse():
            for c in filter(None, r['CAS'].split('|')):
                idx.setdefault(c, r)
        self.assertEqual(idx['108-88-3']['TWA_ppm'], 50.0)       # toluene
        self.assertEqual(idx['1330-20-7']['TWA_ppm'], 100.0)     # xylene all isomers
        self.assertEqual(idx['25551-13-7']['TWA_ppm'], 25.0)     # TMB mixed isomers
        self.assertEqual(idx['78-93-3']['TWA_ppm'], 200.0)       # MEK
        self.assertEqual(idx['107-98-2']['TWA_ppm'], 100.0)      # PGME (multi-line name)


class TestTables(unittest.TestCase):
    def test_every_species_has_properties_and_oel_row(self):
        props = {r['CAS'] for r in rows(V2 / 'data' / 'pure_properties.csv')}
        oel = {r['CAS'] for r in rows(V2 / 'data' / 'oel_assignment.csv')}
        self.assertEqual(props, oel)
        for r in rows(V2 / 'data' / 'pure_properties.csv'):
            self.assertGreater(float(r['Psat25_kPa']), 0)
            self.assertGreater(float(r['Psat35_kPa']), float(r['Psat25_kPa']))
            self.assertTrue(r['Psat_method'])

    def test_surrogate_intervals_feasible(self):
        groups = defaultdict(list)
        for r in rows(V2 / 'data' / 'surrogates.csv'):
            groups[r['UVCB_CAS']].append((float(r['wt_lo']), float(r['wt_hi'])))
        for k, v in groups.items():
            self.assertLessEqual(sum(a for a, _ in v), 100, k)
            self.assertGreaterEqual(sum(b for _, b in v), 100, k)

    def test_tierB_never_weaker_coverage_than_tierA(self):
        for r in rows(V2 / 'data' / 'oel_assignment.csv'):
            if r['TierA_TWA_ppm']:
                self.assertEqual(r['TierA_TWA_ppm'], r['TierB_TWA_ppm'])


class TestComposition(unittest.TestCase):
    def test_central_feasible(self):
        c = run_v2.central([10, 20, 30], [40, 50, 60])
        self.assertAlmostEqual(c.sum(), 100)
        self.assertTrue(((c >= [10, 20, 30]) & (c <= [40, 50, 60])).all())

    def test_range_repair(self):
        lo, hi, note = run_v2.feasible_intervals([60, 50], [70, 60])
        self.assertTrue(note)
        self.assertLessEqual(sum(lo), 100)
        self.assertGreaterEqual(sum(hi), 100)

    def test_expand_conserves_mass(self):
        parts = [{'CAS': '1330-20-7', '함량하한_pct': '40', '함량상한_pct': '60'},
                 {'CAS': '108-88-3', '함량하한_pct': '40', '함량상한_pct': '60'}]
        surr = defaultdict(list)
        for r in rows(V2 / 'data' / 'surrogates.csv'):
            surr[r['UVCB_CAS']].append(r)
        cols, mass, _ = run_v2.expand(parts, surr, 10, random.Random(1))
        self.assertTrue(np.allclose(mass.sum(axis=1), 100))
        self.assertIn(('108-38-3', '1330-20-7'), cols)


class TestFlammability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import run_flammability as rf
        import build_flammability as bf
        cls.rf, cls.bf = rf, bf
        cls.props = {r['CAS']: r for r in rows(V2 / 'data' / 'pure_properties.csv')}
        cls.fl = {r['CAS']: r for r in rows(V2 / 'data' / 'pure_flammability.csv')}
        cls.grids = rf.psat_grid({c: cls.props[c] for c in ('108-88-3', '67-64-1', '112-34-5')})

    def _ratio(self, cas_list):
        fp = np.array([float(self.fl[c]['FP_C']) + 273.15 for c in cas_list])
        return np.array([np.exp(np.log(self.grids[c]) - self.rf.lnp_at(self.grids[c], fp[k]))
                         for k, c in enumerate(cas_list)])

    def test_pure_component_returns_its_own_flash_point(self):
        r = self._ratio(['108-88-3'])
        tf = self.rf.liaw_flash(np.array([[1.0]]), r)[0] - 273.15
        self.assertAlmostEqual(tf, float(self.fl['108-88-3']['FP_C']), delta=1.5)

    def test_ideal_mixture_flash_between_component_values(self):
        cas = ['67-64-1', '112-34-5']          # acetone (-20 C) + 2-(2-butoxyethoxy)ethanol (100 C)
        tf = self.rf.liaw_flash(np.array([[0.5, 0.5]]), self._ratio(cas))[0] - 273.15
        self.assertGreater(tf, float(self.fl['67-64-1']['FP_C']))
        self.assertLess(tf, float(self.fl['112-34-5']['FP_C']))

    def test_crowl_louvar_hexane_close_to_literature(self):
        self.assertAlmostEqual(self.bf.crowl_louvar_pct({'C': 6, 'H': 14}), 1.2, delta=0.15)

    def test_flammability_table_complete(self):
        self.assertEqual(set(self.fl), set(self.props))
        for r in self.fl.values():
            self.assertTrue(r['FP_C'] and r['LFL_pct'] and r['FP_source'], r['CAS'])


class TestCore(unittest.TestCase):
    def test_rank_flip_share(self):
        import pandas as pd
        import analyze_core as ac
        g = pd.Series({'a': 'g', 'b': 'g', 'c': 'g'})
        flips, total = ac.rank_flip_share({'a': 1, 'b': 2, 'c': 3}, {'a': 3, 'b': 2, 'c': 1}, g)
        self.assertEqual((flips, total), (3, 3))
        flips, total = ac.rank_flip_share({'a': 1, 'b': 2, 'c': 3}, {'a': 1, 'b': 2, 'c': 3}, g)
        self.assertEqual((flips, total), (0, 3))

    def test_petroleum_class(self):
        import analyze_flammability as af
        self.assertEqual(af.petroleum_class(20.9), '제1석유류')
        self.assertEqual(af.petroleum_class(21.0), '제2석유류')
        self.assertEqual(af.petroleum_class(70.0), '제3석유류 이상')


if __name__ == '__main__':
    unittest.main()
