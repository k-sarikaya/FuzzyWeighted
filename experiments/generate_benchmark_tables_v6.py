"""Regenerate the two benchmark result tables directly from the current pipeline.

The tables carried over from the earlier submission were produced by a campaign
that predates the correction of the tie-breaking rule, and their reliability
column is stale: the costs are right in all 60 cells, but 24 cells report a
reliability lower than the one the corrected solver attains at the same cost.
This script rebuilds both tables from the solver that is actually described in
the paper, and adds the necessity and possibility reliabilities of each
prescribed design, which the measure family makes available at no extra cost and
which the earlier tables omitted.

    PYTHONPATH=src python3 experiments/generate_benchmark_tables_v6.py
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import ceil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS                      # noqa: E402
from fuzzy_wkon.cpp_backend import solve_at_integer_threshold_cpp  # noqa: E402

ALPHAS = (0.0, 0.3, 0.5, 0.7, 0.9, 1.0)
GAMMAS = (0.0, 0.3, 0.5, 0.7, 1.0)
SPREAD = 0.1
MEASURES = json.loads((ROOT / "results" / "generated" / "measure_family.json").read_text())

TITLES = {
    "pharma_adapted": ("Exact global optimum over the $(\\alpha,\\gamma)$ grid, adapted "
                       "Pharmaceutical Plant catalog", "tab:pharma_stats"),
    "fyffe_adapted":  ("Exact global optimum over the $(\\alpha,\\gamma)$ grid, adapted "
                       "Fyffe 14-subsystem catalog", "tab:fyffe_stats"),
}


def effective_threshold(tfn, alpha, gamma, spread=SPREAD) -> int:
    a, g, d = Fraction(str(alpha)), Fraction(str(gamma)), Fraction(str(spread))
    k1, k2, k3 = (Fraction(str(v)) for v in tfn)
    lower = k1 + a * (k2 - k1)
    upper = k3 - a * (k3 - k2)
    return ceil((upper - g * (upper - lower)) / (1 + d * (1 - a) * (2 * g - 1)))


def measure_cell(catalog: str, alpha: float, gamma: float) -> dict | None:
    for c in MEASURES["catalogs"][catalog]["cells"]:
        if abs(c["alpha"] - alpha) < 1e-9 and abs(c["gamma"] - gamma) < 1e-9:
            return c
    return None


def build(catalog: str) -> tuple[str, list[dict]]:
    benchmark = BENCHMARKS[catalog]
    caption, label = TITLES[catalog]
    rows, body = [], []
    for gamma in GAMMAS:
        for alpha in ALPHAS:
            T = effective_threshold(benchmark["K_tfn"], alpha, gamma)
            out = solve_at_integer_threshold_cpp(benchmark, T, time_limit=900.0)
            mc = measure_cell(catalog, alpha, gamma)
            feasible = out.get("status") == "completed" and out.get("feasible")
            row = {
                "alpha": alpha, "gamma": gamma, "T": T, "feasible": bool(feasible),
                "cost": out.get("best_cost"),
                "R_hurwicz": out.get("best_reliability"),
                "R_necessity": (mc or {}).get("R_necessity"),
                "R_possibility": (mc or {}).get("R_possibility"),
                "nodes": out.get("nodes"), "leaves": out.get("leaves"),
                "seconds": out.get("runtime_seconds"),
                "design": out.get("best_chromosome"),
            }
            rows.append(row)
            if feasible:
                rn = "---" if row["R_necessity"] is None else f"{row['R_necessity']:.3f}"
                rp = "---" if row["R_possibility"] is None else f"{row['R_possibility']:.3f}"
                body.append(f"{alpha:.1f} & {gamma:.1f} & {T} & {row['cost']} & "
                            f"{row['R_hurwicz']:.4f} & {rn} & {rp} & "
                            f"{row['nodes']:,} & {row['seconds']:.2f} \\\\")
            else:
                body.append(f"{alpha:.1f} & {gamma:.1f} & {T} & "
                            f"\\multicolumn{{6}}{{c}}{{$\\dagger$ no feasible design}} \\\\")
        if gamma != GAMMAS[-1]:
            body.append("\\hline")

    tex = (
        "\\begin{table}[!htbp]\n\\centering\n"
        f"\\caption{{{caption}. $T$ is the effective integer threshold of Theorem~\\ref*{{M-thm:reduction}} of the main article; "
        "$R_{\\mathrm{sys}}$ is the reliability of the prescribed design under the Hurwicz rule at that cell, "
        "and $R_N$, $R_\\Pi$ are the necessity and possibility reliabilities of the same design. "
        "$\\dagger$ marks a cell in which the exact solver certifies that no design meets the target.}\n"
        f"\\label{{{label}}}\n\\small\n"
        "\\begin{tabular}{ccrrrrrrr}\n\\hline\n"
        "$\\alpha$ & $\\gamma$ & $T$ & Cost & $R_{\\mathrm{sys}}$ & $R_N$ & $R_\\Pi$ & Nodes & Time (s)\\\\\n"
        "\\hline\n" + "\n".join(body) + "\n\\hline\n\\end{tabular}\n\\end{table}\n"
    )
    return tex, rows


def main() -> None:
    report = {}
    for catalog, filename in (("pharma_adapted", "table_pharma.tex"),
                              ("fyffe_adapted", "table_fyffe.tex")):
        tex, rows = build(catalog)
        (ROOT / "manuscript" / filename).write_text(tex, encoding="utf-8")
        report[catalog] = rows
        feas = [r for r in rows if r["feasible"]]
        print(f"{catalog}: {len(feas)}/{len(rows)} feasible cells written to {filename}")
    dest = ROOT / "results" / "generated" / "benchmark_tables_v6.json"
    dest.write_text(json.dumps(report, indent=2))
    print("written:", dest)


if __name__ == "__main__":
    main()
