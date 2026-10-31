"""Replication table of the main article (v10): the same quantities for the RTS-79 case study and the
two literature-derived instances.

Reads only archived artifacts; it runs no optimization:
  results/generated/benchmark_power_rts.json          30-cell grid, RTS-79
  results/generated/benchmark_tables_v6.json          30-cell grid, Fyffe- and Pharmaceutical-derived
  results/generated/parametric_stability_rts.json     optimal design for every integer threshold, RTS-79
  results/generated/parametric_stability.json         the same for the two derived instances

A design is "at no grid cell" when it is optimal for some integer threshold of the parametric solution
but is the optimum of none of the 30 grid cells.  The script checks that every grid optimum belongs to
the parametric library (Theorem 1), so the difference is well defined.

What the grid misses is also measured on the policy square.  The region on which the effective
threshold is T is, at each alpha < 1, the gamma interval [gamma_T(alpha), gamma_{T-1}(alpha)) of
Corollary 1, so every area below is a one-dimensional integral over alpha (midpoint rule, N_ALPHA
nodes).  On the part of the square whose optimum is at no grid cell the script reports
  * the share of that part on which the design of the nearest grid cell (the cell whose alpha and
    gamma are closest) does not meet the reliability target, or the nearest cell has no design;
  * the largest relative cost excess of the cheapest grid design that does meet the target there;
  * the share of the square on which no grid design meets the target.
Meeting the target is checked with the crisp reliability of the design at the integer threshold
(Theorem 1), with the tolerance 1e-12 of the solver.  These are evaluations of known designs, not
new optimizations.

Writes results/generated/replication_v10.json and manuscript/table_replication.tex.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "generated"
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS          # noqa: E402
from fuzzy_wkon.reliability import static_reliability  # noqa: E402

SPREAD = 0.10
N_ALPHA = 200_000
GRID_ALPHA = (0.0, 0.3, 0.5, 0.7, 0.9, 1.0)
GRID_GAMMA = (0.0, 0.3, 0.5, 0.7, 1.0)


def load(name: str):
    return json.loads((RES / name).read_text(encoding="utf-8"))


def gamma_level(k_tfn, alpha, level):
    """Corollary 1: the gamma at which kappa(alpha, gamma) = level (alpha < 1)."""
    k1, k2, k3 = k_tfn
    lower = k1 + alpha * (k2 - k1)
    upper = k3 - alpha * (k3 - k2)
    return (upper - level * (1.0 - SPREAD * (1.0 - alpha))) / (
        (upper - lower) + 2.0 * level * SPREAD * (1.0 - alpha))


def threshold_area(k_tfn, T, a0=0.0, a1=1.0, g0=0.0, g1=1.0):
    """Area of {(alpha, gamma) in [a0,a1] x [g0,g1] : ceil(kappa) = T}; kappa decreases in gamma."""
    alpha = a0 + (np.arange(N_ALPHA) + 0.5) * (a1 - a0) / N_ALPHA
    lo = np.clip(gamma_level(k_tfn, alpha, T), g0, g1)
    hi = np.clip(gamma_level(k_tfn, alpha, T - 1), g0, g1)
    return float(np.maximum(hi - lo, 0.0).mean() * (a1 - a0))


def bands(nodes):
    """Nearest-node bands of a one-dimensional grid on [0, 1]."""
    cuts = [0.0] + [(nodes[i] + nodes[i + 1]) / 2 for i in range(len(nodes) - 1)] + [1.0]
    return list(zip(cuts[:-1], cuts[1:]))


def meets_target(bench, design, T):
    units = [(bench["choices"][i][j][0], int(bench["choices"][i][j][2]))
             for i, j in enumerate(design) if j != 0]
    return static_reliability(units, T) + 1e-12 >= bench["R_target"]


def design_cost(bench, design):
    return sum(bench["choices"][i][j][1] for i, j in enumerate(design))


def grid_gap(bench, grid_rows, design_key, path):
    """Area measures of what the 30-cell grid misses (see the module docstring)."""
    k_tfn = bench["K_tfn"]
    by_T = {r["T"]: r for r in path}
    assert abs(sum(threshold_area(k_tfn, T) for T in by_T) - 1.0) < 1e-6, "thresholds must tile the square"
    cells = {(r["alpha"], r["gamma"]): r for r in grid_rows}
    grid_designs = {tuple(r[design_key]) for r in grid_rows if r["feasible"]}
    missed = [T for T, r in by_T.items() if r["feasible"] and tuple(r["chromosome"]) not in grid_designs]
    area = {T: threshold_area(k_tfn, T) for T in missed}
    excess, no_grid = [], 0.0
    for T in missed:
        costs = [design_cost(bench, d) for d in grid_designs if meets_target(bench, d, T)]
        if costs:
            excess.append((min(costs) - by_T[T]["cost"]) / by_T[T]["cost"])
        else:
            no_grid += area[T]
    nearest_fails = 0.0
    for a, (a0, a1) in zip(GRID_ALPHA, bands(GRID_ALPHA)):
        for g, (g0, g1) in zip(GRID_GAMMA, bands(GRID_GAMMA)):
            cell = cells[(a, g)]
            for T in missed:
                if not (cell["feasible"] and meets_target(bench, cell[design_key], T)):
                    nearest_fails += threshold_area(k_tfn, T, a0, a1, g0, g1)
    missed_area = sum(area.values())
    return {
        "area_at_no_grid_cell": missed_area,
        "nearest_cell_fails_share": nearest_fails / missed_area,
        "max_cost_excess_cheapest_grid_design": max(excess),
        "area_no_grid_design_meets_target": no_grid,
    }


def summarize(grid_rows, design_key, path, kappa_min, kappa_max, space):
    feasible = [r for r in grid_rows if r["feasible"]]
    grid_designs = {tuple(r[design_key]) for r in feasible}
    par_designs = {tuple(r["chromosome"]) for r in path if r["feasible"]}
    assert grid_designs <= par_designs, "a grid optimum is missing from the parametric library"
    costs = [r["cost"] for r in feasible]
    return {
        "design_space": space,
        "kappa_min": kappa_min,
        "kappa_max": kappa_max,
        "integer_thresholds": len(path),
        "feasible_thresholds": sum(1 for r in path if r["feasible"]),
        "designs_on_square": len(par_designs),
        "designs_at_grid_cells": len(grid_designs),
        "designs_at_no_grid_cell": len(par_designs - grid_designs),
        "grid_cost_min": min(costs),
        "grid_cost_max": max(costs),
        "feasible_cells": len(feasible),
        "cells": len(grid_rows),
    }


def main() -> dict:
    rts_par = load("parametric_stability_rts.json")
    der = load("parametric_stability.json")["catalogs"]
    tables = load("benchmark_tables_v6.json")
    rts_grid = load("benchmark_power_rts.json")
    out = {
        "rts": summarize(rts_grid, "design", rts_par["parametric_path"],
                         rts_par["kappa_min"], rts_par["kappa_max"], "3^10"),
        "fyffe": summarize(tables["fyffe_adapted"], "design", der["fyffe_adapted"]["parametric_path"],
                           der["fyffe_adapted"]["kappa_min"], der["fyffe_adapted"]["kappa_max"], "4^14"),
        "pharma": summarize(tables["pharma_adapted"], "design", der["pharma_adapted"]["parametric_path"],
                            der["pharma_adapted"]["kappa_min"], der["pharma_adapted"]["kappa_max"], "3^10"),
    }
    out["rts"].update(grid_gap(BENCHMARKS["power_rts_adapted"], rts_grid, "design", rts_par["parametric_path"]))
    for key, name in (("fyffe", "fyffe_adapted"), ("pharma", "pharma_adapted")):
        out[key].update(grid_gap(BENCHMARKS[name], tables[name], "design", der[name]["parametric_path"]))
    # consistency with the summary counts stored by the parametric scripts
    assert out["rts"]["designs_on_square"] == rts_par["n_distinct_designs"]
    for key, name in (("fyffe", "fyffe_adapted"), ("pharma", "pharma_adapted")):
        assert out[key]["designs_on_square"] == der[name]["n_distinct_designs"]
        assert out[key]["integer_thresholds"] == len(der[name]["integer_thresholds"])
    (RES / "replication_v10.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    space = {"3^10": r"$3^{10} = 59{,}049$", "4^14": r"$4^{14} = 2.68\times10^{8}$"}
    cols = [out["rts"], out["fyffe"], out["pharma"]]

    def row(label, cells):
        return label + " & " + " & ".join(cells) + r"\\"

    def rng(c):
        return f"$[{c['kappa_min']:.2f},\\,{c['kappa_max']:.2f}]$"

    lines = [
        r"\begin{table}[!htbp]",
        r"\centering",
        r"\caption{The same quantities for the case study and the two literature-derived instances "
        r"($30$-cell grid and complete parametric solution of Theorem~\ref{thm:reduction}). A design ``at no grid cell'' is "
        r"optimal on some region of the policy square but at none of the $30$ grid cells; the last four rows "
        r"describe the region on which such designs are optimal: its area, the share of it on which the design "
        r"of the nearest grid cell (closest $\alpha$ and $\gamma$) does not meet the reliability target, the "
        r"share of the square on which no grid design meets the target, and the largest cost excess over the "
        r"optimum of the cheapest grid design that meets it. Costs are in the "
        r"units of each catalog (M\$ for RTS-79).}",
        r"\label{tab:replication}",
        r"\small",
        r"\begin{tabular}{lccc}",
        r"\hline",
        r" & IEEE RTS-79 & Fyffe-derived & Pharmaceutical-derived\\",
        r"\hline",
        row("Design space", [space[c["design_space"]] for c in cols]),
        row(r"Range of $\kappa(\alpha,\gamma)$", [rng(c) for c in cols]),
        row("Integer thresholds (with a feasible design)",
            [f"{c['integer_thresholds']} ({c['feasible_thresholds']})" for c in cols]),
        row("Distinct optimal designs on the square", [str(c["designs_on_square"]) for c in cols]),
        row("Distinct optimal designs at the grid cells", [str(c["designs_at_grid_cells"]) for c in cols]),
        row("Optimal designs at no grid cell", [str(c["designs_at_no_grid_cell"]) for c in cols]),
        row("Minimum cost on the grid", [f"{c['grid_cost_min']}--{c['grid_cost_max']}" for c in cols]),
        row("Grid cells with a feasible design", [f"{c['feasible_cells']}/{c['cells']}" for c in cols]),
        r"\hline",
        row(r"Area optimal at no grid cell (\% of square)",
            [f"{100 * c['area_at_no_grid_cell']:.1f}" for c in cols]),
        row(r"\quad nearest cell misses target (\% of that area)",
            [f"{100 * c['nearest_cell_fails_share']:.0f}" for c in cols]),
        row(r"\quad no grid design meets target (\% of square)",
            [f"{100 * c['area_no_grid_design_meets_target']:.1f}" for c in cols]),
        row(r"\quad max.\ cost excess of grid design (\%)",
            [f"{100 * c['max_cost_excess_cheapest_grid_design']:.1f}" for c in cols]),
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ]
    (ROOT / "manuscript" / "table_replication.tex").write_text("\n".join(lines), encoding="utf-8")
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
