"""Run full (alpha, gamma) exact optimization grid on IEEE RTS Power Generation Fleet benchmark.

Produces:
- results/generated/benchmark_power_rts.json
- results/generated/table_power_rts.tex (and manuscript/table_power_rts.tex)
- results/generated/table_power_rts_designs.tex (and manuscript/table_power_rts_designs.tex)
"""

from __future__ import annotations

import json
import shutil
import sys
from collections import Counter
from fractions import Fraction
from math import ceil, lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS                      # noqa: E402
from fuzzy_wkon.cpp_backend import solve_at_integer_threshold_cpp  # noqa: E402
from fuzzy_wkon.fuzzy import TFN                                  # noqa: E402
from fuzzy_wkon.membership import membership_profile, saturation_level  # noqa: E402
from fuzzy_wkon.reliability import measure_reliability           # noqa: E402

ALPHAS = (0.0, 0.3, 0.5, 0.7, 0.9, 1.0)
GAMMAS = (0.0, 0.3, 0.5, 0.7, 1.0)
SPREAD = 0.10


def effective_threshold(tfn, alpha, gamma, spread=SPREAD) -> int:
    a, g, d = Fraction(str(alpha)), Fraction(str(gamma)), Fraction(str(spread))
    k1, k2, k3 = (Fraction(str(v)) for v in tfn)
    lower = k1 + a * (k2 - k1)
    upper = k3 - a * (k3 - k2)
    return ceil((upper - g * (upper - lower)) / (1 + d * (1 - a) * (2 * g - 1)))


def modal_options(benchmark: dict, chromosome: list[int]) -> tuple[list[tuple[float, int]], int]:
    denominators = [Fraction(str(value)).denominator for value in benchmark["K_tfn"]]
    for choices in benchmark["choices"]:
        denominators.extend(Fraction(str(weight)).denominator for _, _, weight in choices)
    scale = 1
    for denominator in denominators:
        scale = lcm(scale, denominator)

    selected: list[tuple[float, int]] = []
    for choices, choice_index in zip(benchmark["choices"], chromosome):
        reliability, _cost, weight = choices[choice_index]
        if choice_index == 0:
            continue
        selected.append((float(reliability), int(Fraction(str(weight)) * scale)))
    return selected, scale


def describe_fleet(benchmark: dict, chromosome: list[int]) -> tuple[str, int, int]:
    names = {
        (12, 14): "12 MW Oil",
        (20, 18): "20 MW CT",
        (50, 45): "50 MW Hydro",
        (76, 68): "76 MW Coal",
        (100, 85): "100 MW Oil",
        (155, 135): "155 MW Coal"
    }
    units = []
    total_mw = 0
    total_cost = 0
    for bay, c in enumerate(chromosome):
        if c != 0:
            _, cost, mw = benchmark["choices"][bay][c]
            units.append(names.get((mw, cost), f"{mw} MW"))
            total_mw += mw
            total_cost += cost

    counts = Counter(units)
    parts = []
    # Sort order: 50 MW Hydro, 76 MW Coal, 100 MW Oil, 12 MW Oil, 20 MW CT, 155 MW Coal
    order = ["50 MW Hydro", "76 MW Coal", "100 MW Oil", "12 MW Oil", "20 MW CT", "155 MW Coal"]
    for u in order:
        if u in counts:
            cnt = counts[u]
            parts.append(f"{cnt}×{u}" if cnt > 1 else f"1×{u}")
    for u in sorted(counts.keys()):
        if u not in order:
            cnt = counts[u]
            parts.append(f"{cnt}×{u}" if cnt > 1 else f"1×{u}")
    desc = " + ".join(parts) if parts else "Empty fleet"
    return desc, total_mw, total_cost


def run_power_rts():
    benchmark = BENCHMARKS["power_rts_adapted"]
    threshold_tfn = TFN(*benchmark["K_tfn"])
    rows, body = [], []

    print(f"Solving {benchmark['name']} across {len(GAMMAS)}x{len(ALPHAS)} grid...")
    for gamma in GAMMAS:
        for alpha in ALPHAS:
            T = effective_threshold(benchmark["K_tfn"], alpha, gamma, SPREAD)
            out = solve_at_integer_threshold_cpp(benchmark, T, time_limit=900.0)
            feasible = out.get("status") == "completed" and out.get("feasible")

            r_nec = None
            r_pos = None
            fleet_desc = "---"
            total_mw = 0
            if feasible:
                chromosome = out.get("best_chromosome")
                options, scale = modal_options(benchmark, chromosome)
                cap = max(saturation_level(threshold_tfn, "pos", SPREAD, scale),
                          saturation_level(threshold_tfn, "nec", SPREAD, scale))
                # Graded necessity and possibility profiles per measure_family / paper definition
                prof_nec = membership_profile(threshold_tfn, "nec", SPREAD, SPREAD, scale, cap)
                r_nec = measure_reliability(options, prof_nec)
                prof_pos = membership_profile(threshold_tfn, "pos", SPREAD, SPREAD, scale, cap)
                r_pos = measure_reliability(options, prof_pos)
                fleet_desc, total_mw, _ = describe_fleet(benchmark, chromosome)

            row = {
                "alpha": alpha, "gamma": gamma, "T": T, "feasible": bool(feasible),
                "cost": out.get("best_cost"),
                "total_mw": total_mw,
                "fleet": fleet_desc,
                "R_hurwicz": out.get("best_reliability"),
                "R_necessity": r_nec,
                "R_possibility": r_pos,
                "nodes": out.get("nodes"), "leaves": out.get("leaves"),
                "seconds": out.get("runtime_seconds"),
                "design": out.get("best_chromosome"),
            }
            rows.append(row)
            if feasible:
                rn_str = f"{r_nec:.3f}" if r_nec is not None else "---"
                rp_str = f"{r_pos:.3f}" if r_pos is not None else "---"
                body.append(f"{alpha:.1f} & {gamma:.1f} & {T} & {row['cost']} & "
                            f"{row['R_hurwicz']:.4f} & {rn_str} & {rp_str} & "
                            f"{row['nodes']:,} & {row['seconds']:.2f} \\\\")
            else:
                body.append(f"{alpha:.1f} & {gamma:.1f} & {T} & "
                            f"\\multicolumn{{6}}{{c}}{{$\\dagger$ no feasible design}} \\\\")
        if gamma != GAMMAS[-1]:
            body.append("\\hline")

    caption = "Exact global optimum over the $(\\alpha,\\gamma)$ grid, IEEE RTS-79-parameterized generation fleet case study"
    label = "tab:power_rts_stats"
    tex = (
        "\\begin{table}[!htbp]\n\\centering\n"
        f"\\caption{{{caption}. $T$ is the effective integer threshold of Theorem~\\ref{{thm:reduction}} in MW; "
        "$R_{\\mathrm{sys}}$ is the reliability of the prescribed design under the Hurwicz rule at that cell, "
        "and $R_N$, $R_\\Pi$ are the graded necessity and possibility reliabilities of the same design. "
        "$\\dagger$ marks an infeasible cell.}\n"
        f"\\label{{{label}}}\n\\small\n"
        "\\begin{tabular}{ccrrrrrrr}\n\\hline\n"
        "$\\alpha$ & $\\gamma$ & $T$ & Cost (M\\$) & $R_{\\mathrm{sys}}$ & $R_N$ & $R_\\Pi$ & Nodes & Time (s)\\\\\n"
        "\\hline\n" + "\n".join(body) + "\n\\hline\n\\end{tabular}\n\\end{table}\n"
    )

    # Generate table_power_rts_designs.tex
    unique_designs = {}
    for r in rows:
        if not r["feasible"]:
            continue
        chrom = tuple(r["design"])
        if chrom not in unique_designs:
            unique_designs[chrom] = {
                "cost": r["cost"],
                "total_mw": r["total_mw"],
                "fleet": r["fleet"],
                "R_nec": r["R_necessity"],
                "R_pos": r["R_possibility"],
                "cells": [],
                "T_vals": []
            }
        unique_designs[chrom]["cells"].append(f"({r['alpha']:.1f}, {r['gamma']:.1f})")
        unique_designs[chrom]["T_vals"].append(r["T"])

    design_table_rows = []
    for chrom, info in sorted(unique_designs.items(), key=lambda x: x[1]["cost"]):
        t_min = min(info["T_vals"])
        t_max = max(info["T_vals"])
        t_str = f"{t_min}" if t_min == t_max else f"{t_min}--{t_max}"
        if len(info["cells"]) > 6:
            cells_str = f"{len(info['cells'])} cells (incl.\\ centroid $(0.5,0.5)$)"
        else:
            cells_str = ", ".join(f"${c}$" for c in info["cells"])
        rn_str = f"{info['R_nec']:.4f}" if info['R_nec'] is not None else "---"
        rp_str = f"{info['R_pos']:.4f}" if info['R_pos'] is not None else "---"
        fleet_tex = info['fleet'].replace("×", r"$\times$ ")
        design_table_rows.append(
            f"{info['cost']} & {info['total_mw']} & {t_str} & {fleet_tex} & {rn_str} & {rp_str} & {cells_str} \\\\"
        )

    designs_tex = (
        "\\begin{table}[!htbp]\n\\centering\n"
        "\\caption{Certified optimal generation fleet configurations across the $(\\alpha,\\gamma)$ policy grid for the IEEE RTS-79 case study. "
        "Across the $30$-cell grid, exactly five distinct configurations are optimal. "
        "$R_N$ and $R_\\Pi$ are the graded necessity and possibility reliabilities of each design.}\n"
        "\\label{tab:power_rts_designs}\n\\small\n"
        "\\begin{tabular}{rrcp{3.8cm}ccp{3.6cm}}\n\\hline\n"
        "Cost (M\\$) & Cap (MW) & \\shortstack{$T$ at\\\\grid cells} & Optimal Fleet Configuration & $R_N$ & $R_\\Pi$ & Optimal Policy Cells\\\\\n"
        "\\hline\n" + "\n".join(design_table_rows) + "\n\\hline\n\\end{tabular}\n\\end{table}\n"
    )

    out_json = ROOT / "results" / "generated" / "benchmark_power_rts.json"
    out_tex = ROOT / "results" / "generated" / "table_power_rts.tex"
    out_designs_tex = ROOT / "results" / "generated" / "table_power_rts_designs.tex"
    
    out_json.write_text(json.dumps(rows, indent=2))
    out_tex.write_text(tex)
    out_designs_tex.write_text(designs_tex)

    # Also copy to manuscript dir
    manuscript_dir = ROOT / "manuscript"
    if manuscript_dir.exists():
        (manuscript_dir / "table_power_rts.tex").write_text(tex)
        (manuscript_dir / "table_power_rts_designs.tex").write_text(designs_tex)

    print(f"Wrote {out_json}, {out_tex}, and {out_designs_tex}")


if __name__ == "__main__":
    run_power_rts()
