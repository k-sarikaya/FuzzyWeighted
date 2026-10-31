"""Verification of generalised Proposition 3 under asymmetric proportional spreads (WP-6 / M3).

Verifies that for triangular capacities with asymmetric proportional spreads
w~_i = (w_i(1 - d_l), w_i, w_i(1 + d_r)),
the Hurwicz representative factors as w_i * m(alpha, gamma; d_l, d_r),
where m is strictly independent of component index i.
Tested across 180 random (d_l, d_r) pairs using exact Fraction arithmetic.
"""

from __future__ import annotations

from fractions import Fraction
import random
import unittest


def hurwicz_factor(alpha: Fraction, gamma: Fraction, delta_l: Fraction, delta_r: Fraction) -> Fraction:
    """Return m = 1 - delta_l * (1 - alpha) + gamma * (delta_l + delta_r) * (1 - alpha)."""
    return 1 - delta_l * (1 - alpha) + gamma * (delta_l + delta_r) * (1 - alpha)


class AsymmetricReductionTests(unittest.TestCase):
    def test_factor_independence_180_trials(self) -> None:
        rng = random.Random(20260924)
        trials = 180

        for _ in range(trials):
            # Generate random rational spreads delta_l in (0, 0.4), delta_r in (0, 0.4)
            dl_num = rng.randint(1, 40)
            dr_num = rng.randint(1, 40)
            delta_l = Fraction(dl_num, 100)
            delta_r = Fraction(dr_num, 100)

            # Generate random alpha and gamma in [0, 1]
            alpha = Fraction(rng.randint(0, 100), 100)
            gamma = Fraction(rng.randint(0, 100), 100)

            expected_m = hurwicz_factor(alpha, gamma, delta_l, delta_r)
            self.assertGreater(expected_m, 0, "m must be positive")

            # Test across 10 components with heterogeneous modal capacities
            modal_weights = [Fraction(rng.randint(1, 50)) for _ in range(10)]

            for w_i in modal_weights:
                # TFN endpoints: left = w_i*(1 - delta_l), mode = w_i, right = w_i*(1 + delta_r)
                cut_lower = w_i * (1 - delta_l) + alpha * (w_i - w_i * (1 - delta_l))
                cut_upper = w_i * (1 + delta_r) - alpha * (w_i * (1 + delta_r) - w_i)

                hurwicz_rep = cut_lower + gamma * (cut_upper - cut_lower)
                predicted_rep = w_i * expected_m

                self.assertEqual(
                    hurwicz_rep,
                    predicted_rep,
                    f"Mismatch for w_i={w_i}, delta=({delta_l}, {delta_r}), alpha={alpha}, gamma={gamma}",
                )

            # Also verify summation over an arbitrary operating subset
            subset_weights = modal_weights[:6]
            sum_of_hurwicz = sum(
                (w * (1 - delta_l) + alpha * (w - w * (1 - delta_l)))
                + gamma * (
                    (w * (1 + delta_r) - alpha * (w * (1 + delta_r) - w))
                    - (w * (1 - delta_l) + alpha * (w - w * (1 - delta_l)))
                )
                for w in subset_weights
            )
            sum_of_modals = sum(subset_weights)
            predicted_sum = sum_of_modals * expected_m

            self.assertEqual(sum_of_hurwicz, predicted_sum)


if __name__ == "__main__":
    unittest.main()
