# -*- coding: utf-8 -*-
import unittest
from analyze_ventilated_scenarios import factors, SETTINGS

class TestAnalyticWellMixedRoom(unittest.TestCase):
    def test_zero_source_is_zero(self):
        self.assertEqual(factors(0,1,2,4),(0,0))
        self.assertEqual(factors(2,1,2,0),(0,0))
    def test_ventilation_reduces_integrated_concentration(self):
        low=factors(2,1,0.5,4)[1]
        high=factors(2,1,8,4)[1]
        self.assertGreater(low,high)
    def test_surface_area_increases_index(self):
        small=factors(.1,.2,2,1)[1]
        big=factors(4,.2,2,1)[1]
        self.assertGreater(big,small)
    def test_eight_hour_mean_not_instantaneous_or_saturated(self):
        for _,area,k,ach,t in SETTINGS:
            end,mean=factors(area,k,ach,t)
            self.assertGreater(mean,0)
            self.assertLess(mean,1)
            self.assertLess(end,1)
    def test_no_ventilation_does_not_make_concentration_exceed_saturation(self):
        end,avg=factors(4,2,0,4)
        self.assertGreater(end,avg)
        self.assertLess(avg,1)
    def test_invalid_assumptions_fail_closed(self):
        for params in ((1,1,-1,4),(1,1,1,9),(1,1,1,1,0)):
            with self.assertRaises(ValueError):
                factors(*params)

if __name__=='__main__':
    unittest.main()
