"""The GA baseline must be a fair opponent, not a straw man."""

from __future__ import annotations

import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.heuristics import run_ga
from fuzzy_wkon.solver import solve_exact


class GABaselineTests(unittest.TestCase):
    def test_respects_the_reliability_target(self) -> None:
        """A reported feasible solution must actually meet the target."""

        benchmark = BENCHMARKS["pharma_adapted"]
        result = run_ga(benchmark, 0.5, 0.5, seed=0)
        self.assertTrue(result.feasible)
        self.assertGreaterEqual(result.best_reliability + 1e-12, benchmark["R_target"])

    def test_never_beats_the_certified_optimum(self) -> None:
        """The GA may tie the exact cost but must never report a lower one."""

        benchmark = BENCHMARKS["pharma_adapted"]
        for alpha, gamma in ((0.0, 0.5), (0.5, 0.5), (0.3, 1.0)):
            exact = solve_exact(benchmark, alpha, gamma)
            heuristic = run_ga(benchmark, alpha, gamma, seed=1)
            if not heuristic.feasible:
                continue
            self.assertGreaterEqual(
                heuristic.best_cost, exact["best_cost"],
                msg=f"GA reported a cost below the certified optimum at ({alpha}, {gamma})")

    def test_reports_infeasibility_rather_than_a_bad_design(self) -> None:
        result = run_ga(BENCHMARKS["pharma_adapted"], 0.0, 0.0, seed=0)
        self.assertFalse(result.feasible)
        self.assertIsNone(result.best_cost)


if __name__ == "__main__":
    unittest.main()
