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


if __name__ == '__main__':
    unittest.main()
