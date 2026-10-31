"""The classical readings of Table tab:nesting (main article) on the IEEE RTS-79 fleet.

Each reading is mapped to its cell of the policy square (Proposition 1, Corollary 1), the effective
integer threshold T = ceil(k^alpha(gamma) / m(alpha, gamma)) is computed in exact rational
arithmetic, and the optimal design is read from the archived parametric solution of the fleet
(results/generated/parametric_stability_rts.json) and labelled with the design menu D1-D7
(results/generated/policy_menu_v8.json).  No optimization is run.

Writes results/generated/nesting_rts_example_v11.json and manuscript/table_nesting_rts_supp.tex.

    PYTHONPATH=src python3 experiments/nesting_rts_example_v11.py
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from math import ceil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "generated"
sys.path.insert(0, str(ROOT / "src"))
from fuzzy_wkon.benchmarks import BENCHMARKS  # noqa: E402

SPREAD = F(1, 10)


def k_hurwicz(k_tfn, alpha, gamma):
    k1, k2, k3 = (F(v) for v in k_tfn)
    lo, up = k1 + alpha * (k2 - k1), k3 - alpha * (k3 - k2)
    return up - gamma * (up - lo)


def factor(alpha, gamma):
    return 1 + SPREAD * (1 - alpha) * (2 * gamma - 1)


def threshold(k_tfn, alpha, gamma):
    kap = k_hurwicz(k_tfn, alpha, gamma) / factor(alpha, gamma)
    return kap, ceil(kap)


def readings():
    """(reading, condition, lambda, cell) in the order of Table tab:nesting."""
    half = F(1, 2)
    rows = [("Crisp model (modal values)", r"$W_S \ge k_2$", None, (F(1), half)),
            ("Interval worst case", r"$W_{S,L}^{0} \ge k_U^{0}$", None, (F(0), F(0))),
            ("Interval best case", r"$W_{S,U}^{0} \ge k_L^{0}$", None, (F(0), F(1)))]
    for lam in (F(7, 10), F(9, 10)):
        rows.append(("Necessity", r"$N(\WS\ge\tilde{k}) \ge %s$" % float(lam), lam, (1 - lam, F(0))))
    lam = F(7, 10)
    rows.append(("Possibility", r"$\Pi(\WS\ge\tilde{k}) \ge %s$" % float(lam), lam, (lam, F(1))))
    for lam in (F(3, 10), F(8, 10)):
        cell = (2 * lam, F(1)) if lam <= half else (2 - 2 * lam, F(0))
        rows.append(("Credibility", r"$\mathrm{Cr}(\WS\ge\tilde{k}) \ge %s$" % float(lam), lam, cell))
    return rows


def main() -> dict:
    bench = BENCHMARKS["power_rts_adapted"]
    k_tfn = bench["K_tfn"]
    path = {r["T"]: r for r in json.loads((RES / "parametric_stability_rts.json").read_text())["parametric_path"]}
    menu = json.loads((RES / "policy_menu_v8.json").read_text())["menu"]
    label = {tuple(d["chromosome"]): d for d in menu}
    out = []
    for name, cond, lam, (a, g) in readings():
        kap, T = threshold(k_tfn, a, g)
        opt = path[T]
        d = label[tuple(opt["chromosome"])]
        assert d["T_range"][0] <= T <= d["T_range"][1] and d["cost"] == opt["cost"]
        out.append({"reading": name, "condition": cond, "lambda": None if lam is None else float(lam),
                    "alpha": float(a), "gamma": float(g), "kappa": float(kap), "T": T,
                    "design": d["label"], "cost": opt["cost"], "R_N": d["R_N"]})
    report = {"description": "Readings of Table tab:nesting mapped to cells, thresholds and designs on the IEEE RTS-79 fleet",
              "K_tfn": list(k_tfn), "spread": float(SPREAD), "rows": out}
    (RES / "nesting_rts_example_v11.json").write_text(json.dumps(report, indent=1) + "\n")
    write_table(out, k_tfn)
    return report


def cell_str(a: float, g: float, crisp: bool) -> str:
    return r"$\alpha=1$" if crisp else f"$({a:g},\\,{g:g})$"


def write_table(rows, k_tfn) -> None:
    k1, k2, k3 = k_tfn
    lines = [r"\begin{table}[htbp]", r"\centering",
             r"\caption{The readings of Table~\ref*{M-tab:nesting} of the main article on the IEEE RTS-79 fleet "
             f"($\\tilde{{k}} = ({k1:g}, {k2:g}, {k3:g})$~MW, capacity spread $\\delta = 0.10$). Each reading is mapped to its "
             r"cell; $\kappa$ and $T = \lceil\kappa\rceil$ are the effective continuous and integer thresholds in MW, and the "
             r"design is read from the complete parametric solution (designs D1--D7 of Table~\ref*{M-tab:policy_menu} "
             r"of the main article). No optimization is solved for this table.}",
             r"\label{tab:supp_nesting_rts}", r"\small", r"\begin{tabular}{llcrrlrr}", r"\hline",
             r"Reading & Condition & Cell & $\kappa$ & $T$ & Design & Cost (M\$) & $R_N$\\", r"\hline"]
    for r in rows:
        crisp = r["reading"].startswith("Crisp")
        lines.append(f"{r['reading']} & {r['condition']} & {cell_str(r['alpha'], r['gamma'], crisp)} & "
                     f"{r['kappa']:.2f} & {r['T']} & {r['design']} & {r['cost']} & {r['R_N']:.3f}\\\\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}"]
    (ROOT / "manuscript" / "table_nesting_rts_supp.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    for r in main()["rows"]:
        print(r["reading"], r["alpha"], r["gamma"], round(r["kappa"], 2), r["T"], r["design"], r["cost"])
