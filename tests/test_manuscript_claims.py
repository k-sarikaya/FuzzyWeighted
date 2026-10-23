"""Numbers and statements in the manuscript text checked against the generated artifacts.

Each test reads the LaTeX source and the JSON written by the experiment scripts, so a number that is
edited by hand without re-running the script, or a script result that changes without the text being
updated, makes the test fail.  A second group guards against wordings that were found to be wrong.
"""

import itertools
import json
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS          # noqa: E402
from fuzzy_wkon.reliability import static_reliability  # noqa: E402

MS = ROOT / "manuscript"
RES = ROOT / "results" / "generated"


def tex(name: str) -> str:
    return (MS / name).read_text(encoding="utf-8")


def load(name: str):
    return json.loads((RES / name).read_text(encoding="utf-8"))


def flat(s: str) -> str:
    """Collapse line breaks so that sentences can be matched across source lines."""
    return re.sub(r"\s+", " ", s)


MAIN = flat(tex("manuscript.tex"))
SUPP = flat(tex("supplementary.tex"))
FRONTIER = load("rts_true_frontier_v9.json")
MENU = load("policy_menu_v8.json")
RTARGET = load("rts_rtarget_sensitivity_v9.json")
DELTA = load("delta_sensitivity_v7.json")
BASE = load("baselines_v11.json")
RTS = BENCHMARKS["power_rts_adapted"]


def units(chromosome):
    return [(RTS["choices"][i][c][0], RTS["choices"][i][c][2]) for i, c in enumerate(chromosome) if c]


class FrontierClaims(unittest.TestCase):
    def test_adequate_frontier_table_matches_artifact(self):
        rows = [p for p in FRONTIER["adequate_frontier"] if p["cost"] <= 203]
        table = tex("table_adequate_frontier.tex")
        self.assertEqual([p["cost"] for p in rows], [180, 190, 194, 198, 203])
        for p in rows:
            self.assertIn(f"{p['cost']} & ", table)
            self.assertIn(f"{p['r_n']:.4f}", table)
            self.assertGreaterEqual(p["r_modal"], 0.95 - 1e-12)

    def test_text_numbers_of_the_adequate_frontier(self):
        by_cost = {p["cost"]: p for p in FRONTIER["adequate_frontier"]}
        self.assertIn(f"$R_N = {by_cost[194]['r_n']:.4f}$ for \\$194\\,M", MAIN)
        self.assertIn(f"$R_N = {by_cost[198]['r_n']:.4f}$ for \\$198\\,M", MAIN)
        self.assertIn(f"against ${by_cost[203]['r_n']:.4f}$ for D7", MAIN)
        self.assertIn(f"starts at D6 ($R_N = {by_cost[180]['r_n']:.4f}$)", MAIN)
        self.assertEqual(sorted(FRONTIER["adequate_frontier_not_in_menu"]), [190, 194, 198])
        self.assertIn("three designs that no Hurwicz cell prescribes", MAIN)
        self.assertIn(f"\\$194\\,M and ${by_cost[198]['r_n']:.3f}$ for \\$198\\,M", MAIN)  # budget inversion
        self.assertIn(f"$R_N = {by_cost[194]['r_n']:.3f}$ for \\$194\\,M", MAIN)

    def test_cheapest_adequate_design_costs_180(self):
        self.assertEqual(FRONTIER["cheapest_adequate"]["cost"], 180)
        self.assertIn("the cheapest design with $R_{\\mathrm{modal}} \\ge 0.95$ costs \\$180\\,M", MAIN)
        self.assertIn("With the adequacy constraint} no design cheaper than \\$180\\,M is admissible", MAIN)

    def test_menu_designs_d1_to_d4_have_zero_crisp_adequacy(self):
        modal = {m["label"]: m for m in FRONTIER["menu_modal"]}
        for label, mw in zip(("D1", "D2", "D3", "D4"), (150, 155, 162, 176)):
            self.assertEqual(modal[label]["r_modal"], 0.0)
            self.assertEqual(modal[label]["total_w"], mw)
        self.assertIn("D1--D4 install less than $180$~MW ($150$, $155$, $162$ and $176$~MW)", MAIN)
        for label in ("D5", "D6", "D7"):
            self.assertTrue(modal[label]["adequate"])

    def test_unconstrained_frontier_statements(self):
        pts = {p["cost"]: p for p in FRONTIER["pareto_frontier"]}
        adequate = set(FRONTIER["unconstrained_frontier_adequate_costs"])
        between = [c for c in pts if 158 <= c <= 203]
        self.assertEqual(sorted(c for c in between if c in adequate), [180, 203])
        self.assertIn(f"$R_N = {pts[170]['r_n']:.4f}$", MAIN)
        self.assertAlmostEqual(0.96 ** 2, 0.9216, places=10)
        rns_171 = sorted({round(d["r_n"], 4) for d in FRONTIER["cost_171_designs"]}, reverse=True)
        self.assertEqual(rns_171[0], 0.3603)
        self.assertIn(0.311, {round(r, 3) for r in rns_171})
        self.assertTrue(all(r < pts[170]["r_n"] for r in rns_171))
        self.assertIn("$R_N = 0.3603$", MAIN)
        self.assertIn("$R_N = 0.3110$", MAIN)


class TargetSweepClaims(unittest.TestCase):
    def rows(self, label):
        return {r["R_target"]: r for r in RTARGET if r["case"].startswith(label)}

    def test_sweep_covers_every_target(self):
        for label in ("Optimistic", "Modal", "Necessity"):
            self.assertEqual(sorted(self.rows(label)), [round(0.90 + 0.01 * i, 2) for i in range(10)])

    def test_text_matches_sweep(self):
        modal, opt, nec = self.rows("Modal"), self.rows("Optimistic"), self.rows("Necessity")
        first = {lab: min(t for t, r in rows.items() if r["tolerates_single_loss"])
                 for lab, rows in (("modal", modal), ("opt", opt), ("nec", nec))}
        self.assertEqual(first, {"modal": 0.97, "opt": 0.98, "nec": 0.98})
        r = modal[0.97]
        self.assertIn(f"(\\${r['cost']}\\,M, ${r['total_mw']}$~MW, $\\Rsys = {r['R_sys']:.4f}$)", SUPP)
        self.assertIn(f"\\${opt[0.98]['cost']}\\,M and ${opt[0.98]['total_mw']}$~MW in the optimistic boundary cell", SUPP)
        self.assertIn(f"\\${nec[0.98]['cost']}\\,M and ${nec[0.98]['total_mw']}$~MW in the necessity cell", SUPP)
        self.assertIn(f"optimal at $\\Rtarget=0.97$ in the modal cell (\\${r['cost']}\\,M) and at $0.98$ in the other two "
                      f"(\\${opt[0.98]['cost']}\\,M and \\${nec[0.98]['cost']}\\,M)", MAIN)
        self.assertEqual([modal[t]["cost"] for t in (0.90, 0.92, 0.93, 0.94, 0.95, 0.96)], [167, 167, 172, 172, 180, 180])
        self.assertIn("\\$167\\,M up to $\\Rtarget = 0.92$, \\$172\\,M at $0.93$--$0.94$ and \\$180\\,M at $0.95$--$0.96$", SUPP)
        for rows in (modal, opt, nec):
            r = rows[0.95]
            self.assertFalse(r["tolerates_single_loss"])

    def test_sweep_against_independent_enumeration(self):
        # exhaustive enumeration of all 3^10 fleets, independent of the C++ branch-and-bound
        best = {}
        for chromo in itertools.product(range(3), repeat=10):
            u = units(chromo)
            if not u:
                continue
            cost = sum(RTS["choices"][i][c][1] for i, c in enumerate(chromo) if c)
            for T in (148, 180, 220):
                rel = static_reliability(u, T)
                for t in (0.92, 0.93, 0.97, 0.98):
                    if rel >= t - 1e-12:
                        key = (T, t)
                        if key not in best or (cost, -rel) < best[key]:
                            best[key] = (cost, -rel)
        rows = {(r["T"], r["R_target"]): r for r in RTARGET}
        for key, (cost, neg_rel) in best.items():
            self.assertEqual(rows[key]["cost"], cost, key)
            self.assertAlmostEqual(rows[key]["R_sys"], -neg_rel, places=4)


class SensitivityClaims(unittest.TestCase):
    def test_modal_cell_under_asymmetric_regimes(self):
        cases = {c["label"].split()[0]: c for c in DELTA["capacity_asymmetric_cases"]}
        self.assertEqual((cases["Summer"]["neutral"]["T"], cases["Summer"]["neutral"]["cost"]), (185, 180))
        self.assertEqual((cases["Extreme"]["neutral"]["T"], cases["Extreme"]["neutral"]["cost"]), (185, 180))
        self.assertEqual((cases["Winter"]["neutral"]["T"], cases["Winter"]["neutral"]["cost"]), (176, 158))
        self.assertEqual((cases["Cool"]["neutral"]["T"], cases["Cool"]["neutral"]["cost"]), (176, 158))
        self.assertEqual((cases["Winter"]["optimistic"]["T"], cases["Winter"]["optimistic"]["cost"]), (141, 130))
        self.assertEqual((cases["Winter"]["pessimistic"]["T"], cases["Winter"]["pessimistic"]["cost"]), (209, 203))
        self.assertIn("$T = 185$~MW (\\$180\\,M) for summer and extreme heat and $T = 176$~MW (\\$158\\,M) for winter and cool climate", MAIN)
        self.assertIn("winter, $(0.05,0.15)$, \\$203\\,M at $T=209$~MW", MAIN)


class UnitSpecificClaims(unittest.TestCase):
    def test_cell_0_03_equal_cost_but_d5_infeasible(self):
        cells = load("rts_unit_specific_v7.json")["cells"]
        c = next(r for r in cells if r["alpha"] == 0.0 and r["gamma"] == 0.3)
        self.assertEqual(c["common"]["cost"], c["unit_specific"]["cost"])
        self.assertEqual(c["common"]["cost"], 180)
        k = 198 - 0.3 * (198 - 162)              # k^0(0.3)
        hydro = 0.3 * 50 + 0.7 * 45              # Hurwicz value of (45, 50, 50) at alpha = 0, gamma = 0.3
        self.assertAlmostEqual(k, 187.2)
        self.assertAlmostEqual(4 * hydro, 186.0)
        self.assertLess(4 * hydro, k)
        self.assertIn("four hydro units supply only $4 \\times (0.3 \\cdot 50 + 0.7 \\cdot 45) = 186$~MW", SUPP)


class ImportanceAndCatalogClaims(unittest.TestCase):
    def test_pharma_capacities(self):
        ph = BENCHMARKS["pharma_adapted"]["choices"]
        first = [row[1][2] for row in ph]
        second = [row[2][2] for row in ph]
        self.assertEqual(first, [9, 7, 5, 9, 9, 10, 6, 5, 8, 6])
        self.assertTrue(all(1 <= b - a <= 2 for a, b in zip(first, second)))
        self.assertIn("$w_i = (9, 7, 5, 9, 9, 10, 6, 5, 8, 6)$", SUPP)
        self.assertIn("one or two units larger", SUPP)


class FrontMatterClaims(unittest.TestCase):
    def abstract(self):
        raw = tex("manuscript.tex")
        start = raw.index(r"\abstract{") + len(r"\abstract{")
        return raw[start:raw.index("\n\n\\keywords{", start)].rstrip().removesuffix("}")

    def test_abstract_length_and_numbers(self):
        a = flat(self.abstract())
        words = re.sub(r"\$([^$]*)\$", lambda m: m.group(1).replace("\\", "").replace("{", "").replace(",", ""), a)
        words = words.replace("\\,", " ").replace("\\new{", "").replace("\\begin{abstract}", "")
        self.assertTrue(150 <= len(words.split()) <= 250)      # FODM guide for authors: 150-250 words
        repl = load("replication_v10.json")
        self.assertEqual([(repl[k]["designs_at_no_grid_cell"], repl[k]["designs_on_square"]) for k in ("rts", "fyffe", "pharma")],
                         [(2, 7), (14, 29), (6, 20)])
        shares = [repl[k]["designs_at_no_grid_cell"] / repl[k]["designs_on_square"] for k in ("rts", "fyffe", "pharma")]
        lo, hi = round(100 * min(shares)), round(100 * max(shares))    # 2/7 -> 29, 14/29 -> 48
        self.assertIn(f"a grid over the indices misses {lo}--{hi}\\% of the optimal designs", a)
        rows = [BASE["decision"][k][m] for k in ("rts", "fyffe", "pharma") for m in ("crisp", "w_only", "k_only")]
        self.assertAlmostEqual(max(r["fails_target"] for r in rows), 0.5, places=6)   # "up to half"
        self.assertTrue(all(r["max_cost_excess"] > 0 for r in rows))                 # "overpay"
        self.assertIn("fail the target on up to half of the index range or overpay", a)
        self.assertIn("choose a design instead of fixing the indices in advance", a)
        self.assertIn("the two indices affect the design only through one integer threshold", a)
        self.assertIn("We propose a dual-fuzzy weighted", a)


class ScalingClaims(unittest.TestCase):
    """Runtime bound stated for the small synthetic instances and the 12/12 reduction check."""

    def test_small_instance_runtime_bound(self):
        rows = load("scaling_campaign_v4.json")["rows"]
        small = [r["wall_seconds"] for r in rows if r["n"] in (8, 10)]
        self.assertEqual(len(small), 6)
        self.assertLessEqual(max(small), 0.012)
        self.assertGreater(max(small), 0.005)
        self.assertEqual(MAIN.count("at most $0.012$~s"), 2)


class ReplicationClaims(unittest.TestCase):
    """v10: the replication table and its paragraph are read from the archived grid and parametric artifacts."""

    def test_table_and_text(self):
        sys.path.insert(0, str(ROOT / "experiments"))
        import replication_table_v10 as rep
        rts_par = load("parametric_stability_rts.json")
        der = load("parametric_stability.json")["catalogs"]
        grid = load("benchmark_tables_v6.json")
        got = {
            "rts": rep.summarize(load("benchmark_power_rts.json"), "design", rts_par["parametric_path"],
                                 rts_par["kappa_min"], rts_par["kappa_max"], "3^10"),
            "fyffe": rep.summarize(grid["fyffe_adapted"], "design", der["fyffe_adapted"]["parametric_path"],
                                   der["fyffe_adapted"]["kappa_min"], der["fyffe_adapted"]["kappa_max"], "4^14"),
            "pharma": rep.summarize(grid["pharma_adapted"], "design", der["pharma_adapted"]["parametric_path"],
                                    der["pharma_adapted"]["kappa_min"], der["pharma_adapted"]["kappa_max"], "3^10"),
        }
        archived = load("replication_v10.json")
        for key in got:
            got[key].update(rep.grid_gap(rep.BENCHMARKS[{"rts": "power_rts_adapted", "fyffe": "fyffe_adapted",
                                                         "pharma": "pharma_adapted"}[key]],
                                         load("benchmark_power_rts.json") if key == "rts"
                                         else grid[{"fyffe": "fyffe_adapted", "pharma": "pharma_adapted"}[key]],
                                         "design",
                                         rts_par["parametric_path"] if key == "rts"
                                         else der[{"fyffe": "fyffe_adapted", "pharma": "pharma_adapted"}[key]]["parametric_path"]))
            for name, value in got[key].items():
                if isinstance(value, float):
                    self.assertAlmostEqual(value, archived[key][name], places=9)
                else:
                    self.assertEqual(value, archived[key][name])
        # the integrated RTS-79 area equals the D2 + D6 shares of the design menu, which are lattice estimates
        menu = {d["label"]: d for d in load("policy_menu_v8.json")["menu"]}
        self.assertAlmostEqual(got["rts"]["area_at_no_grid_cell"],
                               menu["D2"]["area_share"] + menu["D6"]["area_share"], delta=2e-3)
        # independent lattice check of the nearest-cell share (argmin distance instead of bands)
        import numpy as np
        b = rep.BENCHMARKS["pharma_adapted"]
        by_T = {r["T"]: r for r in der["pharma_adapted"]["parametric_path"]}
        cells = {(r["alpha"], r["gamma"]): r for r in grid["pharma_adapted"]}
        designs = {tuple(r["design"]) for r in grid["pharma_adapted"] if r["feasible"]}
        x = (np.arange(400) + 0.5) / 400
        k1, k2, k3 = b["K_tfn"]
        missed = fails = 0
        for a in x:
            lo, up = k1 + a * (k2 - k1), k3 - a * (k3 - k2)
            ca = min(rep.GRID_ALPHA, key=lambda v: abs(v - a))
            for g in x:
                T = int(np.ceil((up - g * (up - lo)) / (1 + rep.SPREAD * (1 - a) * (2 * g - 1)) - 1e-12))
                if not by_T[T]["feasible"] or tuple(by_T[T]["chromosome"]) in designs:
                    continue
                missed += 1
                cell = cells[(ca, min(rep.GRID_GAMMA, key=lambda v: abs(v - g)))]
                fails += not (cell["feasible"] and rep.meets_target(b, cell["design"], T))
        self.assertAlmostEqual(missed / x.size ** 2, got["pharma"]["area_at_no_grid_cell"], delta=5e-3)
        self.assertAlmostEqual(fails / missed, got["pharma"]["nearest_cell_fails_share"], delta=1e-2)
        table = flat(tex("table_replication.tex"))
        self.assertIn("Optimal designs at no grid cell & 2 & 14 & 6", table)
        self.assertIn("Integer thresholds (with a feasible design) & 73 (73) & 31 (31) & 22 (20)", table)
        self.assertIn("Minimum cost on the grid & 130--203 & 14--36 & 92--158", table)
        self.assertEqual([got[k]["designs_at_no_grid_cell"] for k in ("rts", "fyffe", "pharma")], [2, 14, 6])
        self.assertIn("$2$ of $7$, $14$ of $29$ and $6$ of $20$ optimal designs are optimal at no grid cell", MAIN)
        self.assertIn("$73$, $31$ and $22$ integer thresholds yield $7$, $29$ and $20$ distinct optimal designs", MAIN)
        self.assertIn("the minimum cost spans $14$--$36$ and $92$--$158$", MAIN)
        pct = lambda f, k, d: f"{100 * archived[k][f]:.{d}f}"
        cover = ", ".join(f"${pct('area_at_no_grid_cell', k, 1)}\\%$" for k in ("rts", "fyffe"))
        self.assertIn(f"their regions cover {cover} and ${pct('area_at_no_grid_cell', 'pharma', 1)}\\%$ of the square", MAIN)
        self.assertIn(f"On ${pct('nearest_cell_fails_share', 'rts', 0)}\\%$, ${pct('nearest_cell_fails_share', 'fyffe', 0)}\\%$ "
                      f"and ${pct('nearest_cell_fails_share', 'pharma', 0)}\\%$ of these regions", MAIN)
        self.assertIn(f"costs up to ${pct('max_cost_excess_cheapest_grid_design', 'rts', 1)}\\%$, "
                      f"${pct('max_cost_excess_cheapest_grid_design', 'fyffe', 1)}\\%$ and "
                      f"${pct('max_cost_excess_cheapest_grid_design', 'pharma', 1)}\\%$ more than the optimum", MAIN)
        self.assertIn(f"On ${pct('area_no_grid_design_meets_target', 'pharma', 1)}\\%$ of the square of the "
                      "Pharmaceutical-derived instance no grid design meets the target", MAIN)
        self.assertEqual([archived[k]["area_no_grid_design_meets_target"] for k in ("rts", "fyffe")], [0.0, 0.0])
        self.assertIn(r"Area optimal at no grid cell (\% of square) & 4.8 & 22.5 & 8.5", table)
        self.assertIn(r"\quad nearest cell misses target (\% of that area) & 53 & 51 & 65", table)
        # the infeasible thresholds of the Pharmaceutical-derived instance are the two largest
        path = der["pharma_adapted"]["parametric_path"]
        infeasible = sorted(r["T"] for r in path if not r["feasible"])
        self.assertEqual(infeasible, sorted(r["T"] for r in path)[-2:])


class BaselineClaims(unittest.TestCase):
    """v11: the model baselines are recomputed and checked against the archive, the tables and the text."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(ROOT / "experiments"))
        import baselines_v11 as bl
        cls.bl = bl

    def test_decision_baselines_recomputed(self):
        for key in ("rts", "fyffe", "pharma"):
            for m in ("crisp", "w_only", "k_only"):
                got = self.bl.decision_baseline(key, m)
                for name, value in got.items():
                    self.assertAlmostEqual(value, BASE["decision"][key][m][name], places=9, msg=(key, m, name))

    def test_crisp_failure_region_is_the_region_above_200_mw(self):
        # independent of the joint integration: the RTS-79 crisp design D5 fails exactly where the dual
        # threshold exceeds 200 MW, an area given by the single-threshold integral of the replication script
        import replication_table_v10 as rep
        k = RTS["K_tfn"]
        lib = self.bl.library("rts")
        area = sum(rep.threshold_area(k, T) for T in lib if T > 200)
        self.assertAlmostEqual(area, BASE["decision"]["rts"]["crisp"]["fails_target"], places=4)
        for m in ("w_only", "k_only"):
            self.assertAlmostEqual(BASE["decision"]["rts"][m]["fails_target"], area, places=4)

    def test_probabilistic_design_independently(self):
        # q_T by a midpoint grid over (D, M) instead of adaptive quadrature; reliabilities by the shared kernel
        import numpy as np
        k1, k2, k3 = RTS["K_tfn"]
        n = 2000
        d = k1 + (np.arange(n) + 0.5) * (k3 - k1) / n
        m = 0.9 + (np.arange(n) + 0.5) * 0.2 / n
        fd = np.where(d <= k2, (d - k1) / ((k3 - k1) * (k2 - k1)), (k3 - d) / ((k3 - k1) * (k3 - k2))) * 2
        fm = np.where(m <= 1.0, (m - 0.9) / (0.2 * 0.1), (1.1 - m) / (0.2 * 0.1)) * 2
        w = np.outer(fd, fm)
        w /= w.sum()
        T = np.ceil(np.outer(d, 1 / m) - 1e-12).astype(int)
        q = {t: float(w[T == t].sum()) for t in np.unique(T)}
        r = BASE["probabilistic"]["rts"]["triangular"]
        design = r["design"]
        rel = {t: static_reliability(units(design), int(t)) for t in q}
        r_prob = sum(q[t] * rel[t] for t in q)
        self.assertAlmostEqual(r_prob, r["R_prob"], places=3)
        self.assertGreaterEqual(r_prob, RTS["R_target"])
        self.assertAlmostEqual(sum(v for t, v in q.items() if t > 200), r["prob_threshold_above_met"], places=3)
        self.assertAlmostEqual(rel[201], 0.99 ** 4 * 0.98, places=12)
        self.assertEqual(rel[213], 0.0)
        self.assertFalse(r["in_parametric_library"])
        self.assertEqual(r["on_square"]["optimal_cost"], 0.0)

    def test_tables_and_text(self):
        table = flat(tex("table_baselines.tex"))
        pct = lambda v: f"{100 * v:.1f}"  # noqa: E731
        for key in ("rts", "fyffe", "pharma"):
            c = BASE["decision"][key]["crisp"]
            self.assertIn(f"{pct(c['fails_target'])} & {pct(c['meets_at_higher_cost'])} & {pct(c['max_cost_excess'])}", table)
        crisp = [BASE["decision"][k]["crisp"] for k in ("rts", "fyffe", "pharma")]
        self.assertIn("The crisp design fails the reliability target on " + ", ".join(
            f"${pct(r['fails_target'])}\\%$" for r in crisp[:2]) + f" and ${pct(crisp[2]['fails_target'])}\\%$", MAIN)
        w = [BASE["decision"][k]["w_only"] for k in ("rts", "fyffe", "pharma")]
        self.assertIn("with fuzzy capacities only the design fails on " + ", ".join(
            f"${pct(r['fails_target'])}\\%$" for r in w[:2]) + f" and ${pct(w[2]['fails_target'])}\\%$", MAIN)
        p = BASE["probabilistic"]["rts"]["triangular"]
        self.assertEqual((p["cost"], p["largest_threshold_met"]), (194, 200))
        self.assertIn(f"at \\$194\\,M ($R_{{\\mathrm{{prob}}}}={p['R_prob']:.4f}$)", MAIN)
        self.assertIn(f"exceeds $200$~MW is ${p['prob_threshold_above_met']:.3f}$", MAIN)
        self.assertIn(f"$R_N={p['necessity']['R_N']:.3f}$, is below the ${p['necessity']['adequate_frontier_R_N_at_budget']:.3f}$", MAIN)
        self.assertIn(f"($R_{{\\mathrm{{modal}}}}={p['necessity']['R_modal']:.3f}$)", MAIN)
        u = BASE["probabilistic"]["rts"]["uniform"]
        self.assertEqual(u["cost"], 198)
        ph = BASE["probabilistic"]["pharma"]["triangular"]
        self.assertIn(f"probabilistic design (cost ${ph['cost']}$) is likewise optimal at no posture "
                      f"and fails the target on ${pct(ph['on_square']['fails_target'])}\\%$", MAIN)


class WordingGuards(unittest.TestCase):
    """Wordings that were found to be wrong must not come back."""

    FORBIDDEN = [
        "coincides precisely", "possibility cell", "Possibility cell", "extensive metaheuristic",
        "bi-objective", "Possibilistic Birnbaum", "possibilistic Birnbaum", "structural importance",
        "structural Birnbaum", "IFRA", "closed-form lifecycle", "closed forms for dynamic",
        "true combinatorial", "disentangl", "x_j=1", "every design cheaper than the crisp",
        "uniquely attributable", "reported as $0.96", "baseline contractual adequacy",
        "weighted-average", "exact analytical formulations", "Aggressive Possibility",
        "changes three times", "both boundaries clipped", "Pillar", "pillar",
    ]

    def test_main_supplement_tables_and_letters(self):
        sources = {"manuscript": MAIN, "supplement": SUPP}
        for path in sorted(MS.glob("table_*.tex")):
            sources[path.name] = flat(path.read_text(encoding="utf-8"))
        for name, text in sources.items():
            for phrase in self.FORBIDDEN:
                self.assertFalse(phrase in text, f"{phrase!r} in {name}")

    def test_american_spelling(self):
        """K4: American spelling in the prose (labels, file names and citation keys excepted)."""
        british = re.compile(r"\b(optimis|minimis|maximis|normalis|factoris|parameteris|characteris|"
                             r"generalis|discretis|hybridis|formalis|organis|stylis)(e|es|ed|er|ing|ation)\b"
                             r"|\bmodell(ed|ing)\b|\banalys(e|ed)\b|\b(neighbour|favour|favours|centre|artefact)\b", re.I)
        protect = re.compile(r"\\(label|ref|eqref|cite[a-z]*|input|includegraphics|texttt|url|href)(\[[^\]]*\])?\{[^}]*\}")
        sources = {"manuscript": MAIN, "supplement": SUPP}
        for path in sorted(MS.glob("table_*.tex")):
            sources[path.name] = flat(path.read_text(encoding="utf-8"))
        for name, text in sources.items():
            found = british.findall(protect.sub("", text))
            self.assertEqual(found, [], f"British spelling in {name}")

    def test_figure_scripts(self):
        for name in ("run_delta_sensitivity_v7.py", "plot_true_frontier_v9.py"):
            text = (ROOT / "experiments" / name).read_text(encoding="utf-8")
            for phrase in ("Aggressive Possibility", "Risk-Neutral Centroid", "R_N=0.311"):
                self.assertFalse(phrase in text, f"{phrase!r} in {name}")


class AuditRecordClaims(unittest.TestCase):
    """Supplementary Section S2 reports the audit counts; they must match the artifact."""

    def test_factorisation_counts_match_artifact(self):
        fac = load("audits_v6.json")["audits"]["factorisation"]
        cells, comparisons = fac["cells"], fac["comparisons"]
        alphas, gammas = len(fac["grid"]["alphas"]), len(fac["grid"]["gammas"])
        self.assertEqual(cells % (alphas * gammas), 0)
        instances = cells // (alphas * gammas)
        claim = (f"${comparisons:,}$".replace(",", "{,}")
                 + f" comparisons across ${instances} \\times ({alphas} \\times {gammas}) = {cells}$ cells")
        self.assertIn(claim, SUPP)
        self.assertEqual(fac["mismatches"], 0)


if __name__ == "__main__":
    unittest.main()
