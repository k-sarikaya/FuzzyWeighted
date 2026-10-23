"""Design menu, budget inversion and gamma intervals of the IEEE RTS-79 fleet (policy_menu_v8.py)."""

import importlib.util
import json
import unittest
from math import ceil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("policy_menu_v8", ROOT / "experiments" / "policy_menu_v8.py")
pm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pm)


class PolicyMenuTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.regions = pm.load_regions()
        cls.grid = json.loads((pm.RESULTS / "benchmark_power_rts.json").read_text(encoding="utf-8"))

    def test_necessity_reliabilities_match_grid_campaign(self):
        # every design that a grid cell prescribes must reproduce that cell's R_N and R_Pi
        for cell in self.grid:
            i = pm.region_of(cell["alpha"], cell["gamma"], self.regions)
            r_n, r_p = pm.measures(self.regions[i]["chromosome"])
            self.assertAlmostEqual(r_n, cell["R_necessity"], places=10)
            self.assertAlmostEqual(r_p, cell["R_possibility"], places=10)

    def test_designs_below_crisp_cost_have_zero_necessity(self):
        for r in self.regions:
            r_n, _ = pm.measures(r["chromosome"])
            if r["cost"] < 180:
                self.assertEqual(r_n, 0.0)
            else:
                self.assertGreater(r_n, 0.5)

    def test_equal_cost_designs_d5_d6(self):
        d5, d6 = self.regions[4], self.regions[5]
        self.assertEqual(d5["cost"], d6["cost"])
        n5, p5 = pm.measures(d5["chromosome"])
        n6, p6 = pm.measures(d6["chromosome"])
        self.assertGreater(n6, n5)   # D6 has the higher necessity reliability
        self.assertGreater(p5, p6)   # D5 has the higher possibility reliability

    def test_gamma_interval_matches_brute_force(self):
        for alpha in (0.0, 0.25, 0.5, 0.75, 0.9):
            for gamma in (0.0, 0.2, 0.45, 0.7, 1.0):
                lo, hi, i = pm.gamma_interval(alpha, gamma, self.regions)
                self.assertLessEqual(lo, gamma + 1e-12)
                self.assertTrue(gamma < hi or hi == 1.0)
                for k in range(1001):
                    g = k / 1000
                    inside = lo <= g < hi or (hi == 1.0 and g == 1.0)
                    same = pm.region_of(alpha, g, self.regions) == i
                    if abs(g - lo) > 1e-9 and abs(g - hi) > 1e-9:
                        self.assertEqual(inside, same, (alpha, gamma, g, lo, hi))


if __name__ == "__main__":
    unittest.main()
