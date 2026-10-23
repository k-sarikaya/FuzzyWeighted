"""Independent exhaustive check: a unit-specific RTS-79 optimum outside the crisp library."""

import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "verify_unit_specific_bruteforce_v8", ROOT / "experiments" / "verify_unit_specific_bruteforce_v8.py")
bf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bf)


class UnitSpecificBruteForceTests(unittest.TestCase):
    def test_hydro_only_optimum_is_outside_crisp_library(self):
        cost, _rel, design = bf.solve(0.5, 0.7, bf.HYDRO_ONLY)
        self.assertEqual(cost, 171)
        self.assertEqual(list(design), [0, 0, 1, 1, 0, 0, 0, 0, 0, 2])   # 2 x 20 MW CT + 1 x 155 MW coal
        path = json.loads((ROOT / "results" / "generated" / "parametric_stability_rts.json").read_text())
        crisp = {tuple(e["chromosome"]) for e in path["parametric_path"] if e.get("chromosome")}
        self.assertNotIn(tuple(design), crisp)

    def test_necessity_cell_matches_compiled_search(self):
        cost, _rel, _design = bf.solve(0.0, 0.0, bf.UNIT)
        self.assertEqual(cost, 203)


if __name__ == "__main__":
    unittest.main()
