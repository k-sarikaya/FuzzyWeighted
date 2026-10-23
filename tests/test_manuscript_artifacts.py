"""Automated consistency test between manuscript text and generated JSON artifacts.

Guarantees that:
1. Table 7 (scaling campaign) numbers (T, nodes, leaves, cost, time) identically match
   results/generated/scaling_campaign_v4.json.
2. Proposition 1 specifies the correct open-closed interval lambda in (0, 1],
   avoiding the lambda = 0 degeneracy.
3. Section 4.11 and Supplementary S2 (formerly Appendix B) accurately report reduction
   verification across all 12 instances with n <= 14.
"""

import json
import re
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT_PATH = REPO_ROOT / "manuscript" / "manuscript.tex"
SCALING_JSON_PATH = REPO_ROOT / "results" / "generated" / "scaling_campaign_v4.json"
POWER_RTS_JSON_PATH = REPO_ROOT / "results" / "generated" / "benchmark_power_rts.json"
SUPPLEMENT_PATH = REPO_ROOT / "manuscript" / "supplementary.tex"


class TestManuscriptArtifactConsistency(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not MANUSCRIPT_PATH.exists():
            raise unittest.SkipTest(f"Manuscript file not found at {MANUSCRIPT_PATH}")
        cls.tex = MANUSCRIPT_PATH.read_text(encoding="utf-8")
        # The reproducibility appendix now lives in Supplementary Material S2.
        cls.supp = SUPPLEMENT_PATH.read_text(encoding="utf-8") if SUPPLEMENT_PATH.exists() else ""
        if SCALING_JSON_PATH.exists():
            cls.scaling_data = json.loads(SCALING_JSON_PATH.read_text(encoding="utf-8"))
        else:
            cls.scaling_data = {"rows": []}

    def test_table_scaling_matches_json(self):
        """Table 7 (tab:scaling) T, cost, leaves, and runtime match scaling_campaign_v4.json."""
        match = re.search(r"\\label\{tab:scaling\}([\s\S]*?)\\end\{tabular\}", self.tex)
        self.assertIsNotNone(match, f"Table tab:scaling not found in {MANUSCRIPT_PATH.name}")

        table_content = match.group(1)
        # Row format: n & Search space $3^n$ & T & Time (s) & Nodes & Leaves & Optimal cost
        # e.g. 16 & $4.30\times10^{7}$ & 78 & 6.20 & $1.49\times10^{7}$ & 29 & 99\\
        pattern = r"(\d+)\s*&\s*\$[^\$]+\$\s*&\s*(\d+)\s*&\s*([\d\.]+)\s*&\s*\$[^\$]+\$\s*&\s*(\d+)\s*&\s*(\d+)"
        rows = re.findall(pattern, table_content)
        self.assertEqual(len(rows), 5, f"Expected 5 rows in Table 7 (n=12,14,16,18,20), found {len(rows)}")

        json_design = {
            r["n"]: r
            for r in self.scaling_data.get("rows", [])
            if r.get("profile") == "design"
        }

        for n_str, T_str, time_str, leaves_str, cost_str in rows:
            n = int(n_str)
            T = int(T_str)
            time_val = float(time_str)
            leaves = int(leaves_str)
            cost = int(cost_str)

            self.assertIn(n, json_design, f"n={n} not found in scaling JSON for profile 'design'")
            jr = json_design[n]

            # 1. Effective threshold T
            self.assertEqual(
                T,
                jr["effective_threshold"],
                f"Mismatch in T for n={n}: TeX={T} vs JSON={jr['effective_threshold']}",
            )
            # 2. Optimal cost
            self.assertEqual(
                cost,
                jr["best_cost"],
                f"Mismatch in Cost for n={n}: TeX={cost} vs JSON={jr['best_cost']}",
            )
            # 3. Evaluated leaves
            self.assertEqual(
                leaves,
                jr["leaves"],
                f"Mismatch in Leaves for n={n}: TeX={leaves} vs JSON={jr['leaves']}",
            )
            # 4. Wall clock time (rounded to 2 decimal places)
            expected_time = round(jr["wall_seconds"], 2)
            self.assertAlmostEqual(
                time_val,
                expected_time,
                places=2,
                msg=f"Mismatch in runtime for n={n}: TeX={time_val} vs JSON={expected_time}",
            )

    def test_proposition_1_domain_lambda(self):
        """Proposition 1 claim must state lambda in (0, 1] to exclude the lambda=0 boundary."""
        prop_block_match = re.search(
            r"\\textbf\{\\stmt\{proposition\}\{prop:chance\}\{Proposition\} \([^\)]+\)\.\}([\s\S]*?\$\\square\$)",
            self.tex,
        )
        self.assertIsNotNone(prop_block_match, "Proposition 1 block not found in the manuscript")
        prop_block = prop_block_match.group(1)

        # Statement must have (0, 1] and not [0, 1]
        claim_match = re.search(r"^([\s\S]*?)\\textit\{Proof\.\}", prop_block)
        self.assertIsNotNone(claim_match, "Statement not found in Proposition 1 block")
        claim_text = claim_match.group(1)
        self.assertIn(r"\lambda \in (0,1]", claim_text)
        self.assertNotIn(r"\lambda \in [0,1]", claim_text)

        # Proof must explicitly explain the restriction
        proof_match = re.search(r"\\textit\{Proof\.\}([\s\S]*?)\$\\square\$", prop_block)
        self.assertIsNotNone(proof_match, "Proof not found in Proposition 1 block")
        proof_text = proof_match.group(1)
        self.assertIn(r"\lambda>0", proof_text)
        self.assertIn(r"\lambda=0", proof_text)

    def test_reduction_check_reporting_consistency(self):
        """Verify that scaling_campaign_v4.json has 12 instances with n<=14, all 12 evaluated, and manuscript agrees."""
        rows_le14 = [r for r in self.scaling_data.get("rows", []) if r.get("n", 999) <= 14]
        self.assertEqual(len(rows_le14), 12, "Expected exactly 12 instances with n <= 14 (sizes 8,10,12,14 x 3 profiles)")

        evaluated_rows = [r for r in rows_le14 if "reduction_check" in r and r["reduction_check"] is not None]
        self.assertEqual(len(evaluated_rows), 12, "Expected exactly 12 evaluated instances with reduction_check")
        for r in evaluated_rows:
            self.assertTrue(r["reduction_check"].get("agrees"), f"Reduction check failed for n={r['n']} {r['profile']}")

        # Check §5.4 wording
        self.assertIn(
            "Across all $12$ instances up to $n=14$, the reduced and\nunreduced formulations agreed on cost, reliability and design.",
            self.tex,
            "Section 5.4 should accurately refer to all 12 instances up to n=14",
        )

        # Check Supplementary S2 (formerly Appendix B) check (i) wording
        self.assertIn(
            "across all $12$\ninstances with $n \\le 14$, all $12$ instances were checked",
            self.supp,
            "Supplementary S2 check (i) should accurately specify all 12 checked instances with n <= 14",
        )

    def test_narrative_consistency_no_strict_barrier_fallacy(self):
        """Ensure outdated claims are absent and new mathematical enhancements are present."""
        self.assertNotIn("Strict-Barrier", self.tex)
        self.assertNotIn("graceful degradation", self.tex)
        self.assertNotIn("over-conservative", self.tex)
        self.assertNotIn("capital waste", self.tex)
        self.assertNotIn("submitted formulation", self.tex)

        # Confirm new narrative anchors and mathematical results are present
        self.assertIn("consider a regional grid operator that must deliver", self.tex)  # motivating example, closing Section 1.1
        self.assertNotIn(r"\subsection{Motivating Example}", self.tex)
        self.assertIn("necessity- and possibility-based chance constraints", self.tex)
        self.assertIn("rate at which capital expenditure trades against the strength of the", self.tex)
        self.assertIn(r"\stmt{corollary}{cor:partition}{Corollary} (Closed-Form Partition Boundaries)", self.tex)
        self.assertIn(r"\gamma_T(\alpha)", self.tex)
        self.assertIn(r"(2-2\lambda,\,0), & \lambda > 1/2.", self.tex)
        self.assertIn(r"\label{tab:nesting}", self.tex)

        # Confirm the audit checks (vi)-(ix), now in Supplementary S2
        self.assertIn("Nine checks are applied", self.supp)
        self.assertIn("reproducing the step-wise sums and $RM=0.60$", self.supp)
        self.assertIn(r"Theorem~\ref*{M-thm:rel_bound} is audited over $20{,}000$ random prefix--completion pairs", self.supp)
        self.assertIn("Arm-U uncertain-measure bound admissibility is audited over $20{,}000$", self.supp)
        self.assertNotIn("\\appendix", self.tex)

        # Confirm Table 1 exact optimization and fuzzy threshold consistency
        self.assertIn(r"$\tilde{k}$ Fuzzy & No & No & No & No & No & No & No & \textbf{Yes}\\", self.tex)
        self.assertIn(r"Design task & Eval. & Eval. & Eval. & \shortstack{Heuristic\\opt.} & Eval. & \shortstack{Heuristic\\opt.} & Eval. & \shortstack{\textbf{Exact}\\\textbf{opt.}}\\", self.tex)

    def test_v7_industrial_case_study_artifacts(self):
        """Verify v7 manuscript incorporates IEEE RTS-79 case study, delta sensitivity, and citations."""
        if MANUSCRIPT_PATH.name in ("manuscript_v8.tex", "manuscript.tex"):
            self.assertIn(r"\input{table_power_rts.tex}", self.tex)
            self.assertIn(r"fig_delta_sensitivity.pdf", self.tex)
            self.assertIn("IEEERTS1979", self.tex)
            self.assertIn("Instances and Protocol", self.tex)
            self.assertIn("The Findings Are Robust to the Spread Model and the Target", self.tex)

            # Confirm JSON artifact exists and has 30 rows
            self.assertTrue(POWER_RTS_JSON_PATH.exists(), "benchmark_power_rts.json missing")
            rts_data = json.loads(POWER_RTS_JSON_PATH.read_text(encoding="utf-8"))
            self.assertEqual(len(rts_data), 30)

    def test_power_rts_artifact_consistency(self):
        """Verify RTS-79 benchmark JSON numbers: 30 cells, 5 unique designs, key cell values."""
        self.assertTrue(POWER_RTS_JSON_PATH.exists(), "benchmark_power_rts.json missing")
        rts_data = json.loads(POWER_RTS_JSON_PATH.read_text(encoding="utf-8"))
        self.assertEqual(len(rts_data), 30)

        # Unique designs count across 30 cells
        unique_chroms = {tuple(r["design"]) for r in rts_data if r.get("feasible")}
        self.assertEqual(len(unique_chroms), 5, f"Expected exactly 5 unique designs in 30-cell RTS grid, found {len(unique_chroms)}")

        cell_map = {(r["alpha"], r["gamma"]): r for r in rts_data}

        # (0.0, 1.0) Optimistic
        c_opt = cell_map.get((0.0, 1.0))
        self.assertIsNotNone(c_opt)
        self.assertEqual(c_opt["cost"], 130)
        self.assertEqual(c_opt["T"], 148)
        self.assertEqual(c_opt["total_mw"], 150)
        self.assertAlmostEqual(c_opt["R_necessity"], 0.000, delta=1e-3)
        self.assertAlmostEqual(c_opt["R_possibility"], 0.0864, delta=1e-3)

        # (0.5, 0.5) Neutral
        c_neut = cell_map.get((0.5, 0.5))
        self.assertIsNotNone(c_neut)
        self.assertEqual(c_neut["cost"], 180)
        self.assertEqual(c_neut["T"], 180)
        self.assertEqual(c_neut["total_mw"], 200)
        self.assertAlmostEqual(c_neut["R_necessity"], 0.5056, delta=1e-3)
        self.assertAlmostEqual(c_neut["R_possibility"], 0.9641, delta=1e-3)

        # (0.0, 0.0) Conservative
        c_pess = cell_map.get((0.0, 0.0))
        self.assertIsNotNone(c_pess)
        self.assertEqual(c_pess["cost"], 203)
        self.assertEqual(c_pess["T"], 220)
        self.assertEqual(c_pess["total_mw"], 226)
        self.assertAlmostEqual(c_pess["R_necessity"], 0.9509, delta=1e-3)
        self.assertAlmostEqual(c_pess["R_possibility"], 0.9782, delta=1e-3)

    def test_delta_sensitivity_artifact_consistency(self):
        """Verify delta_sensitivity_v7.json numbers."""
        delta_path = REPO_ROOT / "results" / "generated" / "delta_sensitivity_v7.json"
        if not delta_path.exists():
            self.skipTest("delta_sensitivity_v7.json not yet generated")
        d_data = json.loads(delta_path.read_text(encoding="utf-8"))
        sym = {item["delta"]: item for item in d_data["capacity_symmetric_sweep"]}
        self.assertEqual(len(sym), 10)
        self.assertEqual(sym[0.02]["pessimistic"]["cost"], 180)
        self.assertEqual(sym[0.02]["optimistic"]["cost"], 149)
        self.assertEqual(sym[0.10]["pessimistic"]["cost"], 203)
        self.assertEqual(sym[0.10]["optimistic"]["cost"], 130)
        self.assertEqual(sym[0.20]["pessimistic"]["cost"], 265)
        self.assertEqual(sym[0.20]["optimistic"]["cost"], 127)

        # Invariance of neutral centroid across all delta
        for d_val, item in sym.items():
            self.assertEqual(item["neutral"]["cost"], 180, f"Neutral cost not 180 at delta={d_val}")

        asym = {item["label"]: item for item in d_data["capacity_asymmetric_cases"]}
        self.assertEqual(asym["Summer Heatwave (0.15, 0.05)"]["pessimistic"]["cost"], 231)
        self.assertEqual(asym["Extreme Heatwave (0.20, 0.10)"]["pessimistic"]["cost"], 265)
        self.assertEqual(asym["Winter Favorable (0.05, 0.15)"]["pessimistic"]["cost"], 203)
        self.assertEqual(asym["Cool Climate (0.10, 0.20)"]["pessimistic"]["cost"], 203)

    def test_parametric_stability_rts_consistency(self):
        """Verify parametric stability for RTS has 73 thresholds and 7 distinct designs."""
        rts_path = REPO_ROOT / "results" / "generated" / "parametric_stability_rts.json"
        if not rts_path.exists():
            self.skipTest("parametric_stability_rts.json not yet generated")
        rts_param = json.loads(rts_path.read_text(encoding="utf-8"))
        self.assertEqual(len(rts_param["integer_thresholds"]), 73)
        self.assertEqual(rts_param["n_distinct_designs"], 7)
        self.assertEqual(rts_param["integer_thresholds"][0], 148)
        self.assertEqual(rts_param["integer_thresholds"][-1], 220)


if __name__ == "__main__":
    unittest.main()

