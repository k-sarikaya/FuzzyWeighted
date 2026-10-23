"""Unit tests for the uncertainty-theoretic oracle (Wei et al. 2025).

Verifies:
1. Exact reproduction of the published numerical example from Wei et al. (2025).
2. The root-bound value 0.8488 on the Pharmaceutical benchmark that certifies
   its categorical infeasibility under Arm U.
"""

import unittest
from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.uncertain_oracle import (
    evaluate_design_uncertain,
    optimistic_bound_uncertain,
    verify_example,
)


class TestUncertainOracle(unittest.TestCase):
    def test_reproduces_wei_et_al_example(self):
        """Verify published numerical sums and uncertain measure from Wei et al. (2025)."""
        res = verify_example()
        self.assertAlmostEqual(res["published_RM_traverse"], 0.6, places=4)
        self.assertAlmostEqual(res["computed_RM"], 0.6, places=5)
        self.assertAlmostEqual(res["traverse_RM_at_step_0.1"], 0.6, places=4)
        self.assertEqual(res["computed_sums"][0.9], 8.3)

    def test_pharma_root_bound_certifies_infeasibility(self):
        """Assert root bound on Pharmaceutical catalog is ~0.8488, below target 0.85."""
        benchmark = BENCHMARKS["pharma_adapted"]
        threshold = float(benchmark["K_tfn"][1])  # 40.0
        target = benchmark["R_target"]  # 0.85

        root_bound = optimistic_bound_uncertain(benchmark, [], threshold, 0.1, 0.1)

        # Manuscript §5.6 cites 0.8488
        self.assertAlmostEqual(root_bound, 0.8488, places=4)
        self.assertLess(
            root_bound,
            target,
            f"Root bound {root_bound:.4f} must be strictly below target {target} to prove infeasibility at root",
        )

    def test_single_component_evaluation(self):
        """Evaluating an all-maximum design on pharma should match the root bound."""
        benchmark = BENCHMARKS["pharma_adapted"]
        threshold = float(benchmark["K_tfn"][1])
        # Choice that maximizes reliability for each subsystem:
        best_design = []
        for choices in benchmark["choices"]:
            best_j = max(range(len(choices)), key=lambda j: choices[j][0])
            best_design.append(best_j)

        _, rm = evaluate_design_uncertain(benchmark, best_design, threshold, 0.1, 0.1)
        self.assertAlmostEqual(rm, 0.8488, places=4)


if __name__ == "__main__":
    unittest.main()
