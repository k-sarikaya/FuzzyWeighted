"""Unit tests for the IEEE RTS-79 Power Generation Fleet industrial benchmark.

Verifies:
1. Structural integrity and physical parameters (10 generation bays, MW capacities, RTS-79 failure rates).
2. Certified optimality across risk stances using the exact solver (C++ backend).
"""

import unittest
from fuzzy_wkon.benchmarks import BENCHMARKS
try:
    from fuzzy_wkon.cpp_backend import BackendUnavailable, ensure_binary, solve_exact_cpp
    try:
        ensure_binary()
        CPP_AVAILABLE = True
    except BackendUnavailable:
        CPP_AVAILABLE = False
except Exception:
    CPP_AVAILABLE = False


class TestPowerRTSBenchmark(unittest.TestCase):
    def setUp(self):
        self.benchmark = BENCHMARKS["power_rts_adapted"]

    def test_benchmark_structure(self):
        """Verify N=10, K_tfn=(162, 180, 198) MW, R_target=0.95 and valid choices."""
        self.assertEqual(self.benchmark["N"], 10)
        self.assertEqual(self.benchmark["K_tfn"], (162.0, 180.0, 198.0))
        self.assertEqual(self.benchmark["R_target"], 0.95)
        self.assertEqual(len(self.benchmark["choices"]), 10)

        for bay_idx, choices in enumerate(self.benchmark["choices"]):
            self.assertTrue(len(choices) >= 2, f"Bay {bay_idx} has fewer than 2 choices")
            # Choice 0 is uncommitted bay (0, 0, 0)
            self.assertEqual(choices[0][0], 0.0)
            self.assertEqual(choices[0][1], 0)
            for j, (p, c, w) in enumerate(choices):
                if j == 0:
                    continue
                self.assertGreaterEqual(p, 0.85, f"Bay {bay_idx} option {j} has unrealistic low reliability {p}")
                self.assertLessEqual(p, 1.0, f"Bay {bay_idx} option {j} has reliability > 1.0")
                self.assertGreater(c, 0, f"Bay {bay_idx} option {j} has non-positive cost")
                self.assertGreater(w, 0, f"Bay {bay_idx} option {j} has non-positive capacity")

    def test_exact_solutions_across_stances(self):
        """Verify certified optimal solutions for neutral, pessimistic, and optimistic stances."""
        if not CPP_AVAILABLE:
            self.skipTest("C++ backend not available")

        # Neutral stance (0.5, 0.5): T = 180 MW, cost = 180
        res_neutral = solve_exact_cpp(self.benchmark, 0.5, 0.5)
        self.assertTrue(res_neutral["feasible"])
        self.assertEqual(res_neutral["best_cost"], 180)
        self.assertGreaterEqual(res_neutral["best_reliability"], 0.95)

        # Pessimistic stance (0.0, 0.0): T = 220 MW, cost = 203
        res_pessimistic = solve_exact_cpp(self.benchmark, 0.0, 0.0)
        self.assertTrue(res_pessimistic["feasible"])
        self.assertEqual(res_pessimistic["best_cost"], 203)
        self.assertGreaterEqual(res_pessimistic["best_reliability"], 0.95)

        # Optimistic stance (0.0, 1.0): T = 148 MW, cost = 130
        res_optimistic = solve_exact_cpp(self.benchmark, 0.0, 1.0)
        self.assertTrue(res_optimistic["feasible"])
        self.assertEqual(res_optimistic["best_cost"], 130)
        self.assertGreaterEqual(res_optimistic["best_reliability"], 0.95)

    def test_graded_measures_definition(self):
        """Verify graded R_N and R_Pi match measure_family single-definition principle."""
        from fractions import Fraction
        from math import lcm
        from fuzzy_wkon.fuzzy import TFN
        from fuzzy_wkon.membership import membership_profile, saturation_level
        from fuzzy_wkon.reliability import measure_reliability

        tfn = TFN(*self.benchmark["K_tfn"])
        spread = 0.10

        chrom_180 = [0, 0, 2, 2, 1, 1, 0, 0, 0, 0]  # $180M design
        chrom_130 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 1]  # $130M design
        chrom_203 = [0, 0, 0, 2, 1, 1, 0, 1, 0, 0]  # $203M design

        def get_options_and_cap(chrom):
            denominators = [Fraction(str(v)).denominator for v in self.benchmark["K_tfn"]]
            for choices in self.benchmark["choices"]:
                denominators.extend(Fraction(str(w)).denominator for _, _, w in choices)
            scale = 1
            for d in denominators:
                scale = lcm(scale, d)
            options = []
            for choices, ci in zip(self.benchmark["choices"], chrom):
                if ci == 0:
                    continue
                r, _, w = choices[ci]
                options.append((float(r), int(Fraction(str(w)) * scale)))
            cap = max(saturation_level(tfn, "pos", spread, scale),
                      saturation_level(tfn, "nec", spread, scale))
            return options, scale, cap

        # 1. 180 design: R_N = 0.5056 +- 1e-3, R_Pi = 0.9641 +- 1e-3
        opts_180, scale_180, cap_180 = get_options_and_cap(chrom_180)
        prof_nec_180 = membership_profile(tfn, "nec", spread, spread, scale_180, cap_180)
        prof_pos_180 = membership_profile(tfn, "pos", spread, spread, scale_180, cap_180)
        rn_180 = measure_reliability(opts_180, prof_nec_180)
        rp_180 = measure_reliability(opts_180, prof_pos_180)
        self.assertAlmostEqual(rn_180, 0.5056, delta=1e-3)
        self.assertAlmostEqual(rp_180, 0.9641, delta=1e-3)

        # 2. 130 design: R_N = 0.0000, R_Pi = 0.0864 +- 1e-3
        opts_130, scale_130, cap_130 = get_options_and_cap(chrom_130)
        prof_nec_130 = membership_profile(tfn, "nec", spread, spread, scale_130, cap_130)
        prof_pos_130 = membership_profile(tfn, "pos", spread, spread, scale_130, cap_130)
        rn_130 = measure_reliability(opts_130, prof_nec_130)
        rp_130 = measure_reliability(opts_130, prof_pos_130)
        self.assertAlmostEqual(rn_130, 0.0000, delta=1e-4)
        self.assertAlmostEqual(rp_130, 0.0864, delta=1e-3)

        # 3. 203 design: R_N = 0.9509 +- 1e-3, R_Pi = 0.9782 +- 1e-3
        opts_203, scale_203, cap_203 = get_options_and_cap(chrom_203)
        prof_nec_203 = membership_profile(tfn, "nec", spread, spread, scale_203, cap_203)
        prof_pos_203 = membership_profile(tfn, "pos", spread, spread, scale_203, cap_203)
        rn_203 = measure_reliability(opts_203, prof_nec_203)
        rp_203 = measure_reliability(opts_203, prof_pos_203)
        self.assertAlmostEqual(rn_203, 0.9509, delta=1e-3)
        self.assertAlmostEqual(rp_203, 0.9782, delta=1e-3)

    def test_unit_spreads_homogeneous_reproduces_common_solution(self):
        """Assert that homogeneous unit_spreads bitwise matches common left/right_spread."""
        if not CPP_AVAILABLE:
            self.skipTest("C++ backend not available")
        homogeneous_us = [[(0.0, 0.0)] + [(0.10, 0.10)] * (len(bay) - 1) for bay in self.benchmark["choices"]]
        res_common = solve_exact_cpp(self.benchmark, 0.0, 0.0, left_spread=0.10, right_spread=0.10)
        res_unit = solve_exact_cpp(self.benchmark, 0.0, 0.0, unit_spreads=homogeneous_us)
        self.assertEqual(res_unit["best_cost"], res_common["best_cost"])
        self.assertEqual(res_unit["best_chromosome"], res_common["best_chromosome"])
        self.assertAlmostEqual(res_unit["best_reliability"], res_common["best_reliability"], places=7)

    def _unit_spreads(self, hydro):
        spreads = {12: (0.05, 0.02), 20: (0.15, 0.05), 50: hydro, 76: (0.05, 0.02),
                   100: (0.05, 0.02), 155: (0.05, 0.02)}
        return [[(0.0, 0.0)] + [spreads[w] for _p, _c, w in bay[1:]] for bay in self.benchmark["choices"]]

    def test_unit_specific_necessity_cells_with_rts79_hydro(self):
        """Hydro (0.10, 0) from RTS-79 Table 6: (0,0) keeps the $203M fleet; (0.3,0) saves $23M."""
        if not CPP_AVAILABLE:
            self.skipTest("C++ backend not available")
        us = self._unit_spreads((0.10, 0.0))
        res = solve_exact_cpp(self.benchmark, 0.0, 0.0, unit_spreads=us, time_limit=60.0)
        self.assertEqual(res["status"], "completed")
        self.assertEqual(res["best_cost"], 203)
        res = solve_exact_cpp(self.benchmark, 0.3, 0.0, unit_spreads=us, time_limit=60.0)
        self.assertEqual(res["status"], "completed")
        self.assertEqual(res["best_cost"], 180)
        self.assertEqual(res["best_chromosome"], [0, 0, 0, 0, 0, 1, 0, 0, 0, 2])  # 1 hydro + 1 x 155 MW coal

    def test_necessity_cell_hinges_on_hydro_downside(self):
        """At (0,0) four hydro units (200 MW) clear 198 MW only if hydro has no downside."""
        if not CPP_AVAILABLE:
            self.skipTest("C++ backend not available")
        costs = {}
        for dl in (0.0, 0.02):
            res = solve_exact_cpp(self.benchmark, 0.0, 0.0, unit_spreads=self._unit_spreads((dl, 0.0)),
                                  time_limit=60.0)
            self.assertEqual(res["status"], "completed")
            costs[dl] = res["best_cost"]
        self.assertEqual(costs[0.0], 180)
        self.assertEqual(costs[0.02], 203)

if __name__ == "__main__":
    unittest.main()

