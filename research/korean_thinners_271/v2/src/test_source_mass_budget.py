# -*- coding: utf-8 -*-
import unittest
from analyze_source_mass_budget import mass_required_kg,truncate_to_stock_kg

class TestFiniteMassBudget(unittest.TestCase):
    def test_zero_source_and_empty_time(self):
        self.assertEqual(mass_required_kg(0,2,1,1,2),0)
        self.assertEqual(mass_required_kg(2,0,1,1,2),0)
    def test_monotonic_mass_loss(self):
        a=mass_required_kg(1,3,1,1,2)
        b=mass_required_kg(4,3,1,1,2)
        self.assertGreater(b,a)
    def test_small_stock_shorter_use(self):
        t,expected,depleted=truncate_to_stock_kg(4,1,3,1,1,2)
        self.assertTrue(depleted)
        self.assertLess(t,4)
        self.assertAlmostEqual(
            mass_required_kg(t,3,1,1,2),1,places=8)
    def test_sufficient_stock_is_no_depletion(self):
        t,expected,depleted=truncate_to_stock_kg(4,100,3,1,1,2)
        self.assertEqual(t,4)
        self.assertFalse(depleted)
    def test_negative_stock_raises(self):
        with self.assertRaises(ValueError):
            truncate_to_stock_kg(1,-1,1,1,1,2)

if __name__=='__main__': unittest.main()
