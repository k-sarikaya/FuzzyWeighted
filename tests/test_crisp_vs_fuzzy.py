"""Crisp versus dual-fuzzy reliability of the same designs (crisp_vs_fuzzy_v8.py)."""

import importlib.util
import json
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("crisp_vs_fuzzy_v8", ROOT / "experiments" / "crisp_vs_fuzzy_v8.py")
cf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cf)

import sys  # noqa: E402
sys.path.insert(0, str(ROOT / "src"))
from fuzzy_wkon.reliability import static_reliability  # noqa: E402


class CrispVsFuzzyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        param = json.loads((cf.RESULTS / "parametric_stability_rts.json").read_text(encoding="utf-8"))
        cls.regions = [r for r in param["stability_regions"] if r["feasible"]]
        cls.grid = json.loads((cf.RESULTS / "benchmark_power_rts.json").read_text(encoding="utf-8"))

    def test_model_columns_match_grid_campaign(self):
        # "both" at (0,0) and (0,1), and "crisp" at alpha = 1, are the Hurwicz reliabilities of the campaign
        for cell in self.grid:
            if (cell["alpha"], cell["gamma"]) not in ((0.0, 0.0), (0.0, 1.0), (1.0, 0.5)):
                continue
            u = cf.units(cell["design"])
            name = {(0.0, 0.0): "nec_both", (0.0, 1.0): "pos_both", (1.0, 0.5): "crisp"}[(cell["alpha"], cell["gamma"])]
            f, k = cf.CONDITIONS[name]
            self.assertAlmostEqual(float(cf.reliability(u, f, k)), cell["R_hurwicz"], places=10)

    def test_enumeration_matches_capped_convolution(self):
        # independent check of every entry against the solver's convolution at a scaled integer threshold
        scale = 10
        for r in self.regions:
            u = cf.units(r["chromosome"])
            for name, (f, k) in cf.CONDITIONS.items():
                options = [(float(p), int(w * f * scale)) for w, p in u]
                threshold = int(-(-k * scale // 1))
                self.assertTrue(all(F(o[1]) == w * f * scale for o, (w, _) in zip(options, u)))
                self.assertAlmostEqual(float(cf.reliability(u, f, k)), static_reliability(options, threshold), places=12)

    def test_one_sided_fuzzification_misses_the_effect(self):
        d5 = cf.units(self.regions[4]["chromosome"])
        crisp = cf.reliability(d5, *cf.CONDITIONS["crisp"])
        self.assertEqual(cf.reliability(d5, *cf.CONDITIONS["nec_w_only"]), crisp)
        self.assertEqual(cf.reliability(d5, *cf.CONDITIONS["nec_k_only"]), crisp)
        self.assertEqual(cf.reliability(d5, *cf.CONDITIONS["nec_both"]), 0)


if __name__ == "__main__":
    unittest.main()
