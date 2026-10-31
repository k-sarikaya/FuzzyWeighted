"""Selecting the policy parameters on the IEEE RTS-79 fleet: design menu, budget inversion, gamma intervals.

The complete parametric solution (parametric_stability_rts.json) contains a finite library of
designs. For each design this script reports its cost, its graded necessity and possibility
reliabilities R_N and R_Pi (same definitions and code path as run_power_rts_campaign_v7.py), the share
of the policy square on which it is optimal, and the grid cells that prescribe it. It also provides

  * budget_inversion(B): the design with the largest R_N among those costing at most B, and
  * gamma_interval(alpha, gamma): the interval of gamma, at fixed alpha, on which the optimal design
    is the same as at (alpha, gamma), from the closed-form boundaries of Corollary 2.

Writes results/generated/policy_menu_v8.json and manuscript/table_policy_menu.tex.

    python experiments/policy_menu_v8.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from math import ceil
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS                               # noqa: E402
from fuzzy_wkon.fuzzy import TFN                                           # noqa: E402
from fuzzy_wkon.membership import membership_profile, saturation_level    # noqa: E402
from fuzzy_wkon.reliability import measure_reliability                    # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "campaign", ROOT / "experiments" / "run_power_rts_campaign_v7.py")
campaign = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(campaign)

RESULTS = ROOT / "results" / "generated"
SPREAD = 0.10
B = BENCHMARKS["power_rts_adapted"]
K1, K2, K3 = B["K_tfn"]
BUDGETS = (150, 180, 203)


def kappa(alpha, gamma):
    lower = K1 + alpha * (K2 - K1)
    upper = K3 - alpha * (K3 - K2)
    return (upper - gamma * (upper - lower)) / (1.0 + SPREAD * (1.0 - alpha) * (2.0 * gamma - 1.0))


def gamma_boundary(alpha: float, level: float) -> float:
    """Corollary 2: the gamma at which kappa(alpha, gamma) = level."""
    lower = K1 + alpha * (K2 - K1)
    upper = K3 - alpha * (K3 - K2)
    return (upper - level * (1.0 - SPREAD * (1.0 - alpha))) / (
        (upper - lower) + 2.0 * SPREAD * level * (1.0 - alpha))


def load_regions() -> list[dict]:
    param = json.loads((RESULTS / "parametric_stability_rts.json").read_text(encoding="utf-8"))
    return [r for r in param["stability_regions"] if r["feasible"]]


def region_of(alpha: float, gamma: float, regions: list[dict]) -> int:
    T = ceil(kappa(alpha, gamma) - 1e-12)
    return next(i for i, r in enumerate(regions) if r["T_min"] <= T <= r["T_max"])


def measures(chromosome: list[int]) -> tuple[float, float]:
    options, scale = campaign.modal_options(B, chromosome)
    k = TFN(*B["K_tfn"])
    cap = max(saturation_level(k, "pos", SPREAD, scale), saturation_level(k, "nec", SPREAD, scale))
    r_n = measure_reliability(options, membership_profile(k, "nec", SPREAD, SPREAD, scale, cap))
    r_p = measure_reliability(options, membership_profile(k, "pos", SPREAD, SPREAD, scale, cap))
    return r_n, r_p


def gamma_interval(alpha: float, gamma: float, regions: list[dict]) -> tuple[float, float, int]:
    """Gamma interval at fixed alpha < 1 on which the optimal design equals that of (alpha, gamma).

    kappa is decreasing in gamma for alpha < 1, and region i holds kappa in (T_min - 1, T_max], so
    the interval is [gamma_T(T_max), gamma_T(T_min - 1)), clipped to [0, 1].
    """
    i = region_of(alpha, gamma, regions)
    r = regions[i]
    lo = max(0.0, gamma_boundary(alpha, r["T_max"]))
    hi = min(1.0, gamma_boundary(alpha, r["T_min"] - 1))
    return lo, hi, i


def budget_inversion(budget: float, menu: list[dict]) -> dict:
    affordable = [d for d in menu if d["cost"] <= budget]
    best = max(affordable, key=lambda d: (d["R_N"], -d["cost"]))
    return {"budget": budget, "design": best["label"], "cost": best["cost"], "R_N": best["R_N"]}


def main() -> None:
    regions = load_regions()
    grid = json.loads((RESULTS / "benchmark_power_rts.json").read_text(encoding="utf-8"))
    n = 801
    A, G = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    T = np.ceil(kappa(A, G) - 1e-12)
    menu = []
    for i, r in enumerate(regions):
        r_n, r_p = measures(r["chromosome"])
        fleet, mw, _ = campaign.describe_fleet(B, r["chromosome"])
        area = float(((T >= r["T_min"]) & (T <= r["T_max"])).mean())
        cells = [(c["alpha"], c["gamma"]) for c in grid
                 if c.get("feasible") and region_of(c["alpha"], c["gamma"], regions) == i]
        menu.append({"label": f"D{i + 1}", "fleet": fleet, "capacity_mw": mw, "cost": r["cost"],
                     "T_range": [r["T_min"], r["T_max"]], "R_N": r_n, "R_Pi": r_p,
                     "area_share": area, "grid_cells": len(cells), "chromosome": r["chromosome"]})
    dominated = [d["label"] for d in menu
                 if any(e["cost"] <= d["cost"] and e["R_N"] > d["R_N"] + 1e-12 for e in menu)]
    budgets = [budget_inversion(b, menu) for b in BUDGETS]
    intervals = []
    for a, g in ((0.5, 0.5), (0.0, 0.0), (0.3, 1.0), (0.0, 0.3)):
        lo, hi, i = gamma_interval(a, g, regions)
        intervals.append({"alpha": a, "gamma": g, "design": f"D{i + 1}", "gamma_low": lo, "gamma_high": hi})

    out = {"description": __doc__.split("\n\n")[0], "spread": SPREAD, "menu": menu,
           "dominated_in_cost_and_R_N": dominated, "budget_inversion": budgets,
           "gamma_intervals": intervals}
    (RESULTS / "policy_menu_v8.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    rows = []
    for d in menu:
        fleet_tex = d["fleet"].replace("×", r" $\times$ ")
        rows.append(f"{d['label']} & {fleet_tex} & {d['cost']} & {d['T_range'][0]}--{d['T_range'][1]} & "
                    f"{d['R_N']:.3f} & {d['R_Pi']:.3f} & {100 * d['area_share']:.1f} & {d['grid_cells']} \\\\")
    table = (
        "\\begin{table}[!htbp]\n\\providecommand{\\revtablecolor}{}\\revtablecolor\n\\centering\n"
        "\\caption{Design menu of the IEEE RTS-79 fleet: the seven designs of the complete parametric "
        "solution ($\\delta=0.10$), with the necessity and possibility reliabilities $R_N$, $R_\\Pi$ of each "
        "design, the share of the policy square on which it is optimal and the number of the $30$ grid cells "
        "that prescribe it.}\n"
        "\\label{tab:policy_menu}\n\\footnotesize\n\\begin{tabular}{llrcrrrr}\n\\hline\n"
        "Design & Fleet & Cost (M\\$) & \\shortstack{$T$ range\\\\on square} & $R_N$ & $R_\\Pi$ & Area (\\%) & Grid cells\\\\\n\\hline\n"
        + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n\\end{table}\n")
    (ROOT / "manuscript" / "table_policy_menu.tex").write_text(table, encoding="utf-8")

    for d in menu:
        print(f"{d['label']} {d['fleet']:40s} {d['cost']:4d} R_N={d['R_N']:.4f} R_Pi={d['R_Pi']:.4f} "
              f"area={100 * d['area_share']:.1f}% cells={d['grid_cells']}")
    print("dominated:", dominated)
    print("budgets:", budgets)
    print("intervals:", intervals)


if __name__ == "__main__":
    main()
