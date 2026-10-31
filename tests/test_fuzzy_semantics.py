import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from fuzzy_wkon.fuzzy import TFN, alpha_cut, alpha_derivative, gamma_from_certainty_equivalent, hurwicz_value


class FuzzySemanticsTests(unittest.TestCase):
    def test_alpha_one_collapses_to_mode(self):
        self.assertEqual(alpha_cut(TFN(9.0, 10.0, 11.0), 1.0), (10.0, 10.0))

    def test_symmetric_gamma_half_is_alpha_invariant(self):
        tfn = TFN(9.0, 10.0, 11.0)
        values = [hurwicz_value(tfn, alpha, 0.5) for alpha in (0.0, 0.3, 0.7, 1.0)]
        self.assertTrue(all(abs(value - 10.0) < 1e-12 for value in values))

    def test_asymmetric_derivative(self):
        self.assertAlmostEqual(alpha_derivative(TFN(8.0, 10.0, 14.0), 0.5), -1.0)

    def test_certainty_equivalent_gamma(self):
        self.assertAlmostEqual(gamma_from_certainty_equivalent(10.0, 20.0, 15.0), 0.5)


if __name__ == "__main__":
    unittest.main()
