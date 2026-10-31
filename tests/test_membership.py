"""Validation of the possibility/necessity/Hurwicz decision-rule family.

The closed forms in ``fuzzy_wkon.membership`` are never trusted on the strength
of a literature citation alone: they are checked against a brute-force sup-min
grid, and the equivalence between the Hurwicz grid columns and the
possibility/necessity chance constraints is checked exhaustively on randomised
instances.
"""

from __future__ import annotations

import random
import unittest

from fuzzy_wkon.fuzzy import TFN, alpha_cut
from fuzzy_wkon.membership import (
    hurwicz_indicator,
    membership_profile,
    necessity_geq,
    possibility_geq,
    proportional_capacity_tfn,
    saturation_level,
)
from fuzzy_wkon.reliability import (
    measure_reliability,
    static_reliability,
    subset_enumeration_measure_reliability,
)

GRID = 20001


def _triangular_membership(tfn: TFN, x: float) -> float:
    if x < tfn.left or x > tfn.right:
        return 0.0
    if x < tfn.mode:
        return (x - tfn.left) / (tfn.mode - tfn.left) if tfn.mode > tfn.left else 1.0
    if x > tfn.mode:
        return (tfn.right - x) / (tfn.right - tfn.mode) if tfn.right > tfn.mode else 1.0
    return 1.0


def _brute_force_possibility_geq(a: TFN, b: TFN, points: int = GRID) -> float:
    """sup over x >= y of min(mu_a(x), mu_b(y)), evaluated on a dense grid."""

    low = min(a.left, b.left)
    high = max(a.right, b.right)
    grid = [low + (high - low) * index / (points - 1) for index in range(points)]
    mu_a = [_triangular_membership(a, x) for x in grid]
    mu_b = [_triangular_membership(b, y) for y in grid]
    suffix_max = [0.0] * points
    running = 0.0
    for index in range(points - 1, -1, -1):
        running = max(running, mu_a[index])
        suffix_max[index] = running
    return max(min(suffix_max[i], mu_b[i]) for i in range(points))


def _random_tfn(rng: random.Random) -> TFN:
    left = rng.uniform(0.0, 20.0)
    mode = left + rng.uniform(0.0, 8.0)
    right = mode + rng.uniform(0.0, 8.0)
    return TFN(left, mode, right)


class ClosedFormTests(unittest.TestCase):
    def test_possibility_matches_sup_min_grid(self) -> None:
        rng = random.Random(7)
        worst = 0.0
        for _ in range(300):
            a, b = _random_tfn(rng), _random_tfn(rng)
            worst = max(worst, abs(possibility_geq(a, b) - _brute_force_possibility_geq(a, b)))
        self.assertLess(worst, 2e-3, "closed-form possibility disagrees with the sup-min grid")

    def test_necessity_is_the_dual_of_possibility(self) -> None:
        rng = random.Random(13)
        worst = 0.0
        for _ in range(300):
            a, b = _random_tfn(rng), _random_tfn(rng)
            dual = 1.0 - _brute_force_possibility_geq(b, a)
            worst = max(worst, abs(necessity_geq(a, b) - dual))
        self.assertLess(worst, 2e-3, "closed-form necessity disagrees with 1 - Pos(b >= a)")

    def test_necessity_never_exceeds_possibility(self) -> None:
        rng = random.Random(21)
        for _ in range(5000):
            a, b = _random_tfn(rng), _random_tfn(rng)
            self.assertLessEqual(necessity_geq(a, b), possibility_geq(a, b) + 1e-12)


class EquivalenceTests(unittest.TestCase):
    """Proposition: the Hurwicz grid columns *are* chance constraints.

    ``Nec(W >= K) >= lambda`` holds exactly when the Hurwicz indicator accepts
    at ``(alpha, gamma) = (1 - lambda, 0)``, and ``Pos(W >= K) >= lambda``
    holds exactly when it accepts at ``(alpha, gamma) = (lambda, 1)``.
    ``lambda = 0`` is excluded because the chance constraint is then vacuous.
    """

    LEVELS = (0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)

    def _instances(self, seed: int, count: int):
        rng = random.Random(seed)
        for _ in range(count):
            modal = rng.uniform(1.0, 60.0)
            spread = rng.choice([0.05, 0.1, 0.2])
            capacity = proportional_capacity_tfn(modal, spread, spread)
            mode = rng.uniform(5.0, 55.0)
            threshold = TFN(mode - rng.uniform(0.1, 10.0), mode, mode + rng.uniform(0.1, 10.0))
            yield capacity, threshold

    def test_necessity_equals_gamma_zero_column(self) -> None:
        for capacity, threshold in self._instances(11, 4000):
            for level in self.LEVELS:
                self.assertEqual(
                    necessity_geq(capacity, threshold) >= level - 1e-12,
                    hurwicz_indicator(capacity, threshold, 1.0 - level, 0.0) == 1.0,
                )

    def test_possibility_equals_gamma_one_column(self) -> None:
        for capacity, threshold in self._instances(12, 4000):
            for level in self.LEVELS:
                self.assertEqual(
                    possibility_geq(capacity, threshold) >= level - 1e-12,
                    hurwicz_indicator(capacity, threshold, level, 1.0) == 1.0,
                )


class ProfileTests(unittest.TestCase):
    THRESHOLD = TFN(34.0, 40.0, 46.0)

    def test_profiles_are_non_decreasing_and_saturate(self) -> None:
        for measure in ("pos", "nec"):
            cap = saturation_level(self.THRESHOLD, measure, 0.1, 10)
            profile = membership_profile(self.THRESHOLD, measure, 0.1, 0.1, 10, cap)
            self.assertEqual(profile[-1], 1.0)
            self.assertTrue(all(b >= a - 1e-12 for a, b in zip(profile, profile[1:])))

    def test_indicator_profile_reproduces_static_reliability(self) -> None:
        """An indicator profile must return exactly the V2 number."""

        options = [(0.9, 3), (0.8, 2), (0.7, 2)]
        profile = [0.0, 0.0, 0.0, 0.0, 1.0]
        self.assertAlmostEqual(measure_reliability(options, profile),
                               static_reliability(options, 4), places=12)

    def test_graded_kernel_matches_exponential_oracle(self) -> None:
        rng = random.Random(5)
        for _ in range(200):
            size = rng.randint(1, 9)
            options = [(rng.uniform(0.05, 0.99), rng.randint(0, 6)) for _ in range(size)]
            cap = rng.randint(1, 14)
            raw = sorted(rng.uniform(0.0, 1.0) for _ in range(cap))
            profile = raw + [1.0]
            self.assertAlmostEqual(
                measure_reliability(options, profile),
                subset_enumeration_measure_reliability(options, profile),
                places=12,
            )

    def test_sandwich_holds_at_the_reliability_level(self) -> None:
        """R_nec <= R_pos for the same design, by pointwise membership order."""

        options = [(0.88, 7), (0.86, 5), (0.87, 9), (0.85, 10), (0.82, 6)]
        scale = 10
        caps = {m: saturation_level(self.THRESHOLD, m, 0.1, scale) for m in ("pos", "nec")}
        cap = max(caps.values())
        values = {}
        for measure in ("pos", "nec"):
            profile = membership_profile(self.THRESHOLD, measure, 0.1, 0.1, scale, cap)
            values[measure] = measure_reliability(options, profile)
        self.assertLessEqual(values["nec"], values["pos"] + 1e-12)


class CredibilityTests(unittest.TestCase):
    """Verify credibility measure connections (WP-8 / M6)."""

    def test_credibility_modal_equivalence(self) -> None:
        """Cr(A >= B) >= 0.5 <=> A.mode >= B.mode for any continuous TFNs."""
        rng = random.Random(42)
        for _ in range(5000):
            a = _random_tfn(rng)
            b = _random_tfn(rng)
            cr = 0.5 * (possibility_geq(a, b) + necessity_geq(a, b))
            self.assertEqual(cr >= 0.5, a.mode >= b.mode)

    def test_credibility_coincides_with_hurwicz_at_half_for_symmetric_thresholds(self) -> None:
        """Hurwicz indicator at gamma=0.5 matches Cr >= 0.5 for symmetric thresholds."""
        rng = random.Random(123)
        for _ in range(5000):
            modal_w = rng.uniform(10.0, 100.0)
            spread = rng.uniform(0.01, 0.3)
            cap = proportional_capacity_tfn(modal_w, spread, spread)
            # Symmetric threshold
            mode = rng.uniform(20.0, 80.0)
            hw = rng.uniform(1.0, 15.0)
            thresh = TFN(mode - hw, mode, mode + hw)
            alpha = rng.uniform(0.0, 1.0)

            h_val = hurwicz_indicator(cap, thresh, alpha=alpha, gamma=0.5)
            cr_val = 0.5 * (possibility_geq(cap, thresh) + necessity_geq(cap, thresh))
            self.assertEqual(h_val == 1.0, cr_val >= 0.5)

    def test_general_credibility_hurwicz_mapping(self) -> None:
        """Cr(A >= B) >= lambda maps to Hurwicz gamma=0 (for lambda >= 0.5) or gamma=1 (for lambda < 0.5)."""
        rng = random.Random(456)
        for _ in range(5000):
            a = _random_tfn(rng)
            b = _random_tfn(rng)
            lam = rng.uniform(1e-4, 1.0)
            cr = 0.5 * (possibility_geq(a, b) + necessity_geq(a, b))
            cr_pass = (cr >= lam)

            if lam >= 0.5:
                alpha = 2.0 * (1.0 - lam)
                h_pass = (hurwicz_indicator(a, b, alpha=alpha, gamma=0.0) == 1.0)
            else:
                alpha = 2.0 * lam
                h_pass = (hurwicz_indicator(a, b, alpha=alpha, gamma=1.0) == 1.0)

            self.assertEqual(cr_pass, h_pass)


class CredibilityColumnTests(unittest.TestCase):
    """v11 corollary: every credibility chance constraint is a boundary cell of the square.

    Cr(W >= K) >= lambda  <=>  accepted at (2 - 2 lambda, 0)  for lambda in (1/2, 1],
    Cr(W >= K) >= lambda  <=>  accepted at (2 lambda, 1)      for lambda in (0, 1/2],
    with asymmetric capacity spreads and asymmetric thresholds.
    """

    LEVELS = (0.02, 0.1, 0.25, 0.4, 0.5, 0.51, 0.6, 0.75, 0.9, 1.0)

    def test_credibility_maps_to_boundary_columns(self) -> None:
        rng = random.Random(2026)
        checked = 0
        for _ in range(20000):
            modal = rng.uniform(1.0, 60.0)
            capacity = proportional_capacity_tfn(modal, rng.uniform(0.01, 0.3), rng.uniform(0.01, 0.3))
            mode = rng.uniform(5.0, 55.0)
            threshold = TFN(mode - rng.uniform(0.1, 10.0), mode, mode + rng.uniform(0.1, 10.0))
            cr = 0.5 * (possibility_geq(capacity, threshold) + necessity_geq(capacity, threshold))
            for level in self.LEVELS:
                cell = (2.0 - 2.0 * level, 0.0) if level > 0.5 else (2.0 * level, 1.0)
                if abs(cr - level) < 1e-9:
                    continue  # boundary of the level set, decided by rounding only
                self.assertEqual(cr >= level, hurwicz_indicator(capacity, threshold, *cell) == 1.0,
                                 (capacity, threshold, level, cr))
                checked += 1
        self.assertGreater(checked, 190000)


class Corollary3Tests(unittest.TestCase):
    """Verify Corollary 3 behaviour for symmetric vs asymmetric thresholds (WP-4 / M1)."""

    def test_asymmetric_threshold_counterexample(self) -> None:
        """Corollary 3 fails for asymmetric threshold K=(30, 40, 60)."""
        k = TFN(30.0, 40.0, 60.0)
        k_lower_0, k_upper_0 = alpha_cut(k, 0.0)
        kappa_0 = 0.5 * (k_lower_0 + k_upper_0)
        self.assertEqual(kappa_0, 45.0)
        self.assertNotEqual(kappa_0, k.mode)

        k_lower_half, k_upper_half = alpha_cut(k, 0.5)
        kappa_half = 0.5 * (k_lower_half + k_upper_half)
        self.assertEqual(kappa_half, 42.5)
        self.assertNotEqual(kappa_half, k.mode)

        k_lower_1, k_upper_1 = alpha_cut(k, 1.0)
        kappa_1 = 0.5 * (k_lower_1 + k_upper_1)
        self.assertEqual(kappa_1, 40.0)

    def test_symmetric_threshold_holds_for_all_alpha(self) -> None:
        """Corollary 3 holds for symmetric threshold K=(34, 40, 46)."""
        k = TFN(34.0, 40.0, 46.0)
        for alpha in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0):
            k_lower, k_upper = alpha_cut(k, alpha)
            kappa = 0.5 * (k_lower + k_upper)
            self.assertAlmostEqual(kappa, k.mode, places=12)


class Proposition1DegenerateTests(unittest.TestCase):
    """Verify Proposition 1 non-degeneracy condition (WP-5 / M2)."""

    def test_degenerate_crisp_case_reproduced(self) -> None:
        """When both TFNs are crisp at the same mode, (b3-b2)+(a2-a1)=0.
        Nec(A >= B) = 0 while Hurwicz indicator at (0,0) gives 1.0.
        """
        a = TFN(5.0, 5.0, 5.0)
        b = TFN(5.0, 5.0, 5.0)
        nec = necessity_geq(a, b)
        h_ind = hurwicz_indicator(a, b, alpha=0.0, gamma=0.0)
        self.assertEqual(nec, 0.0)
        self.assertEqual(h_ind, 1.0)
        denom = (b.right - b.mode) + (a.mode - a.left)
        self.assertEqual(denom, 0.0)


if __name__ == "__main__":
    unittest.main()

