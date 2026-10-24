# -*- coding: utf-8 -*-
import unittest
import numpy as np

from analyze_regulatory_vapor_gap import status

class TestRegulatoryVaporGap(unittest.TestCase):
    def test_only_below_when_entire_composition_range_below_rule(self):
        self.assertEqual(status(1, 24.99, 25), 'below_rule_threshold')
        self.assertEqual(status(1, 25, 25), 'straddles_threshold')
        self.assertEqual(status(25, 30, 25), 'at_or_above_rule_threshold')
    def test_invalid_range_cannot_be_classified(self):
        self.assertEqual(status(30, 20, 25), 'range_invalid')
        self.assertEqual(status(np.nan, 20, 25), 'range_invalid')
        self.assertEqual(status(0, 120, 25), 'range_invalid')
    def test_cannot_treat_mass_percent_as_ppm(self):
        # The status function only accepts MSDS mass-percent and a SAME-UNIT
        # regulatory mass-percent cutoff; vapor is evaluated separately.
        self.assertEqual(status(20, 30, 85), 'below_rule_threshold')

if __name__ == '__main__':
    unittest.main()
