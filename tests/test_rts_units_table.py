"""Table tab:rts_units of the main article agrees with the RTS-79 catalog in the code."""
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fuzzy_wkon.benchmarks import BENCHMARKS  # noqa: E402

MAIN = (ROOT / "manuscript" / "manuscript.tex").read_text(encoding="utf-8")


class RtsUnitsTable(unittest.TestCase):
    def test_rows_match_catalog(self):
        block = MAIN[MAIN.index(r"\label{tab:rts_units}"):]
        block = block[:block.index(r"\end{tabular}")]
        rows = re.findall(r"^([A-Za-z ]+) & (\d+) & ([\d.]+) & (\d+) & ([\d.]+) & (\d+)--(\d+)\\\\", block, re.M)
        self.assertEqual(len(rows), 6)
        bays = {}
        for i, row in enumerate(BENCHMARKS["power_rts_adapted"]["choices"], start=1):
            for p, c, w in row[1:]:
                bays.setdefault((int(w), round(p, 3), int(c)), []).append(i)
        self.assertEqual(len(bays), 6)
        for _name, w, p, c, ratio, b0, b1 in rows:
            key = (int(w), float(p), int(c))
            self.assertIn(key, bays)
            self.assertEqual((min(bays[key]), max(bays[key])), (int(b0), int(b1)))
            self.assertEqual(sorted(set(bays[key])), list(range(int(b0), int(b1) + 1)))
            self.assertEqual(f"{int(c) / int(w):.2f}", ratio)
        self.assertNotIn("Gas Combustion", MAIN)


if __name__ == "__main__":
    unittest.main()
