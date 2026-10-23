"""The RTS-79 illustration of the nesting table (Supplementary Material, Table tab:supp_nesting_rts).

Thresholds are recomputed here in floating point from the cell of each reading (Proposition 1,
Corollary 1) and compared with the archived JSON and the typeset table; designs are checked
against the parametric library and the design menu.
"""
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "generated"
DATA = json.loads((RES / "nesting_rts_example_v11.json").read_text())
TABLE = (ROOT / "manuscript" / "table_nesting_rts_supp.tex").read_text(encoding="utf-8")
SUPP = (ROOT / "manuscript" / "supplementary.tex").read_text(encoding="utf-8")


def kappa(k, a, g, d=0.10):
    lo, up = k[0] + a * (k[1] - k[0]), k[2] - a * (k[2] - k[1])
    return (up - g * (up - lo)) / (1 + d * (1 - a) * (2 * g - 1))


class NestingExample(unittest.TestCase):
    def test_cells_follow_proposition_1_and_corollary_1(self):
        for r in DATA["rows"]:
            lam = r["lambda"]
            if r["reading"].startswith("Necessity"):
                self.assertEqual((r["alpha"], r["gamma"]), (round(1 - lam, 12), 0.0))
            elif r["reading"].startswith("Possibility"):
                self.assertEqual((r["alpha"], r["gamma"]), (lam, 1.0))
            elif r["reading"].startswith("Credibility"):
                cell = (round(2 * lam, 12), 1.0) if lam <= 0.5 else (round(2 - 2 * lam, 12), 0.0)
                self.assertEqual((r["alpha"], r["gamma"]), cell)

    def test_thresholds_designs_and_table(self):
        k = DATA["K_tfn"]
        lib = {x["T"]: x for x in json.loads((RES / "parametric_stability_rts.json").read_text())["parametric_path"]}
        menu = {d["label"]: d for d in json.loads((RES / "policy_menu_v8.json").read_text())["menu"]}
        for r in DATA["rows"]:
            kap = kappa(k, r["alpha"], r["gamma"])
            self.assertAlmostEqual(kap, r["kappa"], places=9)
            self.assertEqual(math.ceil(kap - 1e-9), r["T"])
            d = menu[r["design"]]
            self.assertTrue(d["T_range"][0] <= r["T"] <= d["T_range"][1])
            self.assertEqual(lib[r["T"]]["cost"], r["cost"])
            self.assertIn(f"{r['kappa']:.2f} & {r['T']} & {r['design']} & {r['cost']} &", TABLE)
        self.assertEqual([r["T"] for r in DATA["rows"]], [180, 220, 148, 208, 216, 170, 167, 203])

    def test_text(self):
        self.assertIn(r"\input{table_nesting_rts_supp.tex}", SUPP)
        self.assertIn("moves the threshold to $208$ or $216$~MW and buys D7 for \\$203\\,M", SUPP)


if __name__ == "__main__":
    unittest.main()
