import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.oracle import exhaustive_optimum, pareto_dp_optimum, vectorized_exhaustive_optimum
from fuzzy_wkon.solver import solve_exact


class SolverTests(unittest.TestCase):
    def test_reviewer_counterexample_is_resolved(self):
        result = solve_exact(BENCHMARKS["pharma_adapted"], 0.3, 1.0)
        self.assertTrue(result["feasible"])
        self.assertEqual(result["best_cost"], 103)
        self.assertGreaterEqual(result["best_reliability"], 0.8836184302 - 1e-12)

    def test_exact_solver_matches_oracle_on_reduced_instance(self):
        source = BENCHMARKS["pharma_adapted"]
        benchmark = {
            "name": "reduced oracle fixture",
            "N": 4,
            "K_tfn": (12.0, 14.0, 16.0),
            "R_target": 0.60,
            "choices": source["choices"][:4],
        }
        for alpha, gamma in ((0.0, 1.0), (0.7, 0.5), (1.0, 1.0)):
            exact = solve_exact(benchmark, alpha, gamma)
            oracle = exhaustive_optimum(benchmark, alpha, gamma)
            self.assertEqual(exact["feasible"], oracle["feasible"])
            self.assertEqual(exact["best_cost"], oracle["cost"])
            self.assertAlmostEqual(exact["best_reliability"], oracle["reliability"], places=12)

    def test_vectorized_oracle_matches_reference_oracle(self):
        source = BENCHMARKS["pharma_adapted"]
        benchmark = {
            "name": "reduced vectorized oracle fixture",
            "N": 4,
            "K_tfn": (12.0, 14.0, 16.0),
            "R_target": 0.60,
            "choices": source["choices"][:4],
        }
        reference = exhaustive_optimum(benchmark, 0.3, 1.0)
        vectorized = vectorized_exhaustive_optimum(benchmark, 0.3, 1.0)
        self.assertEqual(reference["cost"], vectorized["cost"])
        self.assertAlmostEqual(reference["reliability"], vectorized["reliability"], places=12)

    def test_pareto_dp_matches_reduced_oracle(self):
        source = BENCHMARKS["pharma_adapted"]
        benchmark = {
            "name": "reduced DP fixture",
            "N": 4,
            "K_tfn": (12.0, 14.0, 16.0),
            "R_target": 0.60,
            "choices": source["choices"][:4],
        }
        reference = exhaustive_optimum(benchmark, 0.3, 1.0)
        dp = pareto_dp_optimum(benchmark, 0.3, 1.0, max_cost=100)
        self.assertEqual(reference["cost"], dp["cost"])
        self.assertAlmostEqual(reference["reliability"], dp["reliability"], places=12)


if __name__ == "__main__":
    unittest.main()
