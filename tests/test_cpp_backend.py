"""The C++ backend must be a faithful port, not merely a fast one.

These tests are skipped when no C++ compiler is present, so the suite still
runs on a machine without a toolchain; they never lower the bar for a machine
that has one.
"""

from __future__ import annotations

import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.solver import solve_exact

try:
    from fuzzy_wkon.cpp_backend import BackendUnavailable, ensure_binary, solve_exact_cpp
    try:
        ensure_binary()
        AVAILABLE = True
        REASON = ""
    except BackendUnavailable as error:  # pragma: no cover - environment dependent
        AVAILABLE, REASON = False, str(error)
except Exception as error:  # pragma: no cover - environment dependent
    AVAILABLE, REASON = False, str(error)


@unittest.skipUnless(AVAILABLE, f"C++ backend unavailable: {REASON}")
class BackendEquivalenceTests(unittest.TestCase):
    CELLS = ((0.0, 0.0), (0.0, 1.0), (0.3, 1.0), (0.5, 0.5), (1.0, 0.5))

    def test_matches_python_solver_on_pharma(self) -> None:
        benchmark = BENCHMARKS["pharma_adapted"]
        for alpha, gamma in self.CELLS:
            expected = solve_exact(benchmark, alpha, gamma)
            actual = solve_exact_cpp(benchmark, alpha, gamma)
            self.assertEqual(actual["status"], "completed")
            self.assertEqual(bool(actual.get("feasible")), bool(expected["feasible"]),
                             msg=f"feasibility differs at ({alpha}, {gamma})")
            if not expected["feasible"]:
                continue
            self.assertEqual(actual["best_cost"], expected["best_cost"],
                             msg=f"cost differs at ({alpha}, {gamma})")
            self.assertAlmostEqual(actual["best_reliability"], expected["best_reliability"],
                                   places=9, msg=f"reliability differs at ({alpha}, {gamma})")
            self.assertEqual(actual["best_chromosome"], expected["best_chromosome"],
                             msg=f"design differs at ({alpha}, {gamma})")

    def test_reviewer_equal_cost_counterexample(self) -> None:
        """The referee's (0.3, 1.0) case must survive the port."""

        actual = solve_exact_cpp(BENCHMARKS["pharma_adapted"], 0.3, 1.0)
        self.assertEqual(actual["best_cost"], 103)
        self.assertAlmostEqual(actual["best_reliability"], 0.8836184303, places=9)

    def test_timeout_is_reported_not_silently_completed(self) -> None:
        from fuzzy_wkon.scaling import synthetic_benchmark

        result = solve_exact_cpp(synthetic_benchmark(24, "stringent"), 0.3, 0.7,
                                 time_limit=0.05)
        self.assertEqual(result["status"], "timeout")
        self.assertNotIn("best_cost", result)


if __name__ == "__main__":
    unittest.main()
