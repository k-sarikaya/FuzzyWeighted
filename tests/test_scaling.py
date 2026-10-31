import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from fuzzy_wkon.oracle import exhaustive_optimum
from fuzzy_wkon.scaling import synthetic_benchmark
from fuzzy_wkon.solver import solve_exact
from fuzzy_wkon.solver import prepare_options


class ScalingFamilyTests(unittest.TestCase):
    def test_generator_is_deterministic(self):
        self.assertEqual(synthetic_benchmark(8, "design"), synthetic_benchmark(8, "design"))

    def test_declared_decimal_tfn_does_not_create_spurious_scale(self):
        benchmark = synthetic_benchmark(32, "stringent")
        _, _, scale = prepare_options(benchmark, 0.3, 0.7)
        self.assertLessEqual(scale, 1000)

    def test_all_profiles_use_three_declared_choices(self):
        for profile in ("screening", "design", "stringent"):
            benchmark = synthetic_benchmark(8, profile)
            self.assertTrue(all(len(row) == 3 for row in benchmark["choices"]))

    def test_small_scaling_instance_matches_exhaustive_oracle(self):
        benchmark = synthetic_benchmark(6, "screening")
        exact = solve_exact(benchmark, 0.3, 0.7)
        oracle = exhaustive_optimum(benchmark, 0.3, 0.7)
        self.assertEqual(exact["feasible"], oracle["feasible"])
        self.assertEqual(exact["best_cost"], oracle["cost"])
        self.assertAlmostEqual(exact["best_reliability"], oracle["reliability"], places=12)


if __name__ == "__main__":
    unittest.main()
