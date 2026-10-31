"""Statements of the supplement that are derived from archived artifacts or from the code."""
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fuzzy_wkon.benchmarks import BENCHMARKS  # noqa: E402
from fuzzy_wkon.uncertain_oracle import _attained, design_weights_and_states  # noqa: E402

SUPP = (ROOT / "manuscript" / "supplementary.tex").read_text(encoding="utf-8")
GEN = ROOT / "results" / "generated"


class ArmUInfeasibilityCause(unittest.TestCase):
    def test_units_credited_at_the_target_level(self):
        bench = BENCHMARKS["pharma_adapted"]
        # per subsystem, the largest unit whose state measure reaches the target level 0.85
        x = []
        for row in bench["choices"]:
            ok = [j for j, (p, _c, w) in enumerate(row) if p >= 0.85 and w > 0]
            x.append(max(ok, key=lambda j: row[j][2]) if ok else 0)
        self.assertEqual([i + 1 for i, c in enumerate(x) if c], [2, 3, 4, 6, 10])
        weights, states = design_weights_and_states(bench, x, 0.1, 0.1)
        credited = _attained(0.85, weights, states)
        self.assertLess(credited, 40.0)
        self.assertEqual(f"{credited:.2f}", "39.99")
        audit = json.loads((GEN / "audits_v6.json").read_text())["audits"]["uncertain_all_designs"]
        self.assertEqual(audit["best_chromosome"], x)
        self.assertIn("units of five subsystems (2, 3, 4, 6 and 10)", SUPP)
        self.assertIn("$39.99$ against the threshold $40$", SUPP)


class ArtifactStatements(unittest.TestCase):
    def test_unit_specific_solver_time(self):
        d = json.loads((GEN / "rts_unit_specific_v7.json").read_text())
        times = [c["unit_specific"]["runtime_seconds"] for c in d["cells"]]
        self.assertEqual(len(times), 30)
        self.assertEqual(f"{sum(times):,.0f}", "1,363")
        self.assertEqual(round(max(times)), 234)
        self.assertIn("take $1{,}363$~s of solver time in total, the slowest $234$~s", SUPP)

    def test_published_bisection_value(self):
        v = json.loads((GEN / "uncertainty_baseline_v6.json").read_text())["oracle_verification"]
        self.assertEqual(v["published_RM_bisection"], 0.5977)
        self.assertIn("as $0.5977$ by bisection", SUPP)

    def test_interpreted_scaling_sizes(self):
        d = json.loads((GEN / "scaling_campaign.json").read_text())
        self.assertEqual(d["sizes"], [8, 10, 12, 16, 20, 24, 28, 32])
        solved = {r["n"] for r in d["rows"] if r["status"] == "completed"}
        self.assertEqual(max(solved), 12)
        self.assertIn(r"$n \in \{8, 10, 12, 16, 20, 24, 28, 32\}$ with a $12$~s budget", SUPP)

    def test_reduction_check_record(self):
        d = json.loads((GEN / "scaling_campaign_v4.json").read_text())
        rows = [r for r in d["rows"] if r["n"] <= 14]
        self.assertEqual(sum(1 for r in rows if r["reduction_check"]["agrees"]), 12)
        self.assertIn("(12 instances,", d["reduction_verification"])


if __name__ == "__main__":
    unittest.main()
