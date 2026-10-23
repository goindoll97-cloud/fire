# -*- coding: utf-8 -*-
"""Standard-library regression tests for the source-linked 271-document snapshot."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path
from subprocess import run

ROOT = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]

class DatasetAndCalculationTests(unittest.TestCase):
    def test_validate_snapshot(self):
        p = run([sys.executable, str(ROOT / "src" / "validate_271.py")],
                cwd=REPO, capture_output=True, text=True, check=False)
        self.assertEqual(p.returncode, 0, p.stderr + p.stdout)

    def test_run_raoult_scenarios(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = run([sys.executable, str(ROOT / "src" / "run_raoult_271.py"),
                     "--samples", "10", "--seed", "20261023", "--output", tmp],
                    cwd=REPO, capture_output=True, text=True, check=False)
            self.assertEqual(p.returncode, 0, p.stderr + p.stdout)
            def data(path):
                with (Path(tmp) / path).open(encoding="utf-8-sig", newline="") as f:
                    return list(csv.DictReader(f))
            summary = data("raoult_제품별_요약.csv")
            scenarios = data("raoult_시나리오.csv")
            components = data("raoult_성분별.csv")
            self.assertEqual(len(summary), 6)
            self.assertEqual(len(scenarios), 60)
            self.assertEqual(len(data("raoult_계산제외.csv")), 0)
            for row in scenarios:
                self.assertGreater(float(row["총평형증기압_kPa"]), 0.0)
            from collections import defaultdict
            by_scenario = defaultdict(list)
            for row in components:
                by_scenario[(row["제품ID"], row["시나리오ID"])].append(row)
            self.assertEqual(len(by_scenario), 60)
            for values in by_scenario.values():
                self.assertAlmostEqual(sum(float(r["MSDS질량백분율"]) for r in values), 100., places=3)
                self.assertAlmostEqual(sum(float(r["액상몰분율"]) for r in values), 1., places=6)
                self.assertAlmostEqual(sum(float(r["평형기상몰분율"]) for r in values), 1., places=6)

if __name__ == "__main__":
    unittest.main()
