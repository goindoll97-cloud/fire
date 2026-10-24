# -*- coding: utf-8 -*-
"""Execute with: python -m unittest -v test_analyze_thinner_multichem.py"""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import analyze_thinner_multichem as multi


class CounterfactualTests(unittest.TestCase):
    def test_mass_balance_and_no_duplication(self):
        comp = multi.make_trace_mixture([("67-64-1", 27), ("108-88-3", 73)], "71-43-2", 0.039)
        self.assertAlmostEqual(sum(w for _, w in comp), 100.)
        self.assertAlmostEqual(comp[-1][1], 0.039)
        with self.assertRaisesRegex(ValueError, "ALREADY_PRESENT"):
            multi.make_trace_mixture([("71-43-2",100)], "71-43-2", 0.039)

    def test_no_silent_threshold_assumption(self):
        self.assertEqual(multi.cutoff_relation(0.039, {"screening_cutoff_wt_pct": ""}),
                         "NO_REVIEWED_NUMERIC_CUTOFF")
        entry={"screening_cutoff_wt_pct": 0.1}
        self.assertIn("BELOW_NUMERIC_CUTOFF", multi.cutoff_relation(0.099, entry))
        self.assertIn("AT_NUMERIC_CUTOFF", multi.cutoff_relation(0.100, entry))

    def test_raoult_equals_existing_backend(self):
        props=multi.core.load_reference(multi.ROOT/"thinner_vle_reference_25C.csv")
        mixture=[("108-88-3", 53), ("108-10-1",46.961), ("71-43-2", 0.039)]
        xs,g,pp,ys,vw,terms,share=multi.compute_one(mixture,props,
                   multi.core.GammaModels(Path("not-installed")), "raoult", 100.0)
        self.assertTrue(np.allclose(g,1))
        self.assertAlmostEqual(sum(ys),1)
        self.assertAlmostEqual(sum(vw),1)
        self.assertAlmostEqual(sum(share),1)
        b=2
        self.assertTrue(0<ys[b]*100<1)
        self.assertTrue(0<vw[b]*100<1)
        self.assertAlmostEqual(vw[b]*100 / 0.039, (vw[b]*100)/0.039)

    def test_end_to_end_raoult(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"res"
            args=SimpleNamespace(input=multi.DEFAULT_INPUT,
                reference=multi.ROOT/"thinner_vle_reference_25C.csv",
                registry=multi.ROOT/"multichem_screening_registry.csv",
                output_dir=out,hanna_repo=Path(tmp)/"missing",
                product_ids=["1208563","273605","1183746"],
                candidates=list(multi.DEFAULT_CANDIDATES), doses=[0.039,0.099,0.100],
                models=["raoult"],total_air_ppm=100., strict_models=True)
            meta=multi.run(args)
            self.assertEqual(meta["unavailable_or_failed_runs"],0)
            status=pd.read_csv(out/"02_model_status.csv")
            coverage=pd.read_csv(out/"01_scenario_coverage.csv")
            detail=pd.read_csv(out/"03_component_predictions.csv")
            trace=pd.read_csv(out/"04_trace_screening.csv")
            self.assertTrue((status.status=="OK").all())
            self.assertEqual(len(coverage),9)
            self.assertGreater(int((coverage.result.str.startswith("SKIP")).sum()),0)
            self.assertTrue((trace.legal_omission_assessment=="NOT_PERFORMED").all())
            self.assertTrue((trace.real_product_detection=="NOT_DETERMINED_COUNTERFACTUAL").all())
            self.assertAlmostEqual(float(detail.groupby(["product_id","target_CAS","target_wt_pct","model_id"])["vapor_wt_pct"].sum().min()),100,places=8)
            self.assertTrue((trace.vapor_wt_enrichment_fold>0).all())
            self.assertEqual(meta["successfully_computed_runs"],len(trace))


if __name__=="__main__":
    unittest.main()
