import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from fuzzy_wkon.reliability import static_reliability, subset_enumeration_reliability


class ReliabilityTests(unittest.TestCase):
    def test_capped_convolution_matches_independent_oracle(self):
        rng = random.Random(42)
        for _ in range(20):
            options = [(rng.random(), rng.randrange(0, 8)) for _ in range(7)]
            threshold = rng.randrange(1, 15)
            self.assertAlmostEqual(
                static_reliability(options, threshold),
                subset_enumeration_reliability(options, threshold),
                places=12,
            )


if __name__ == "__main__":
    unittest.main()
