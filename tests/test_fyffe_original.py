import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from fuzzy_wkon.fyffe_original import load_table_i, parallel_reliability, solve_original_fyffe


class OriginalFyffeTests(unittest.TestCase):
    def test_table_i_has_the_published_structure(self):
        table = load_table_i()
        self.assertEqual(len(table), 14)
        self.assertEqual(sum(len(unit) for unit in table), 48)
        self.assertEqual((table[0][2].reliability, table[0][2].cost, table[0][2].weight), (0.91, 2, 2))
        self.assertEqual((table[13][3].reliability, table[13][3].cost, table[13][3].weight), (0.99, 6, 9))

    def test_parallel_reliability(self):
        self.assertAlmostEqual(parallel_reliability(0.9, 1), 0.9)
        self.assertAlmostEqual(parallel_reliability(0.9, 2), 0.99)

    def test_reproduces_published_table_iii_solution(self):
        result = solve_original_fyffe(130, 170)
        self.assertEqual(result["best_cost"], 119)
        self.assertEqual(result["best_weight"], 170)
        self.assertAlmostEqual(result["best_reliability"], 0.9700, places=4)


if __name__ == "__main__":
    unittest.main()
