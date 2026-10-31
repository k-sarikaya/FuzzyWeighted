"""The parametric-partition figure: closed-form boundaries and grid consistency (RTS-79)."""

import importlib.util
import unittest
from math import ceil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "plot_parametric_partition_v7", ROOT / "experiments" / "plot_parametric_partition_v7.py")
fig = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fig)


class ParametricPartitionTests(unittest.TestCase):
    def setUp(self):
        self.param, self.grid, self.regions = fig.load()

    def test_boundaries_and_grid_cells_agree(self):
        hit = fig.check(self.param, self.grid, self.regions)   # raises on any mismatch
        self.assertEqual(len(self.param["integer_thresholds"]), 73)
        self.assertEqual(len(self.regions), 7)
        self.assertEqual(sorted(set(hit)), [0, 2, 3, 4, 6])        # D2 and D6 on no grid cell

    def test_designs_missed_by_the_grid(self):
        # D2 on the possibility column for alpha in (1/11, 17/67]; D6 on the necessity
        # column for alpha in [27/77, 9/19)
        def region(alpha, gamma):
            T = ceil(fig.kappa(alpha, gamma) - 1e-12)
            return next(i for i, r in enumerate(self.regions) if r["T_min"] <= T <= r["T_max"])
        self.assertEqual(region(0.10, 1.0), 1)
        self.assertEqual(region(0.25, 1.0), 1)
        self.assertEqual(region(0.36, 0.0), 5)
        self.assertEqual(region(0.47, 0.0), 5)
        self.assertEqual(self.regions[5]["cost"], self.regions[4]["cost"])
        self.assertLess(self.regions[5]["cost"], self.regions[6]["cost"])


if __name__ == "__main__":
    unittest.main()
