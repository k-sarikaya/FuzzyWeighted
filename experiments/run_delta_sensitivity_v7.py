"""Run capacity tolerance and demand forecast sensitivity analysis on IEEE RTS fleet.

Produces:
- results/generated/delta_sensitivity_v7.json
- figures/fig_delta_sensitivity.pdf
- figures/fig_delta_sensitivity.png
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import ceil, lcm
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS                      # noqa: E402
from fuzzy_wkon.cpp_backend import solve_exact_cpp               # noqa: E402
from fuzzy_wkon.fuzzy import TFN                                  # noqa: E402
from fuzzy_wkon.membership import membership_profile, saturation_level  # noqa: E402
from fuzzy_wkon.reliability import measure_reliability           # noqa: E402

DELTAS = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
ETAS = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]

ASYMMETRIC_CASES = [
    {"label": "Summer Heatwave (0.15, 0.05)", "delta_l": 0.15, "delta_r": 0.05},
    {"label": "Extreme Heatwave (0.20, 0.10)", "delta_l": 0.20, "delta_r": 0.10},
    {"label": "Winter Favorable (0.05, 0.15)", "delta_l": 0.05, "delta_r": 0.15},
    {"label": "Cool Climate (0.10, 0.20)", "delta_l": 0.10, "delta_r": 0.20},
]


def describe_fleet(benchmark: dict, chromosome: list[int]) -> str:
    names = {
        (12, 14): "12 MW Oil",
        (20, 18): "20 MW CT",
        (50, 45): "50 MW Hydro",
        (76, 68): "76 MW Coal",
        (100, 85): "100 MW Oil",
        (155, 135): "155 MW Coal"
    }
    from collections import Counter
    units = []
    for bay, c in enumerate(chromosome):
        if c != 0:
            _, cost, mw = benchmark["choices"][bay][c]
            units.append(names.get((mw, cost), f"{mw} MW"))
    counts = Counter(units)
    order = ["50 MW Hydro", "76 MW Coal", "100 MW Oil", "12 MW Oil", "20 MW CT", "155 MW Coal"]
    parts = []
    for u in order:
        if u in counts:
            cnt = counts[u]
            parts.append(f"{cnt}×{u}" if cnt > 1 else f"1×{u}")
    for u in sorted(counts.keys()):
        if u not in order:
            cnt = counts[u]
            parts.append(f"{cnt}×{u}" if cnt > 1 else f"1×{u}")
    return " + ".join(parts) if parts else "Empty"


def compute_measures(benchmark: dict, chromosome: list[int],
                     spread_l: float, spread_r: float) -> tuple[float, float]:
    """Computes graded necessity and possibility reliabilities under the actual spread."""
    threshold_tfn = TFN(*benchmark["K_tfn"])
    # Modal options and scale
    denominators = [Fraction(str(v)).denominator for v in benchmark["K_tfn"]]
    for choices in benchmark["choices"]:
        denominators.extend(Fraction(str(w)).denominator for _, _, w in choices)
    scale = 1
    for d in denominators:
        scale = lcm(scale, d)

    options = []
    for choices, ci in zip(benchmark["choices"], chromosome):
        if ci == 0:
            continue
        r, _, w = choices[ci]
        options.append((float(r), int(Fraction(str(w)) * scale)))

    cap = max(saturation_level(threshold_tfn, "pos", spread_r, scale),
              saturation_level(threshold_tfn, "nec", spread_l, scale))
    prof_nec = membership_profile(threshold_tfn, "nec", spread_l, spread_r, scale, cap)
    r_nec = measure_reliability(options, prof_nec)
    prof_pos = membership_profile(threshold_tfn, "pos", spread_l, spread_r, scale, cap)
    r_pos = measure_reliability(options, prof_pos)
    return r_nec, r_pos


def run_delta_sensitivity():
    base_benchmark = BENCHMARKS["power_rts_adapted"]
    # Demand is strictly fixed for capacity derating sensitivity
    fixed_k_tfn = base_benchmark["K_tfn"]  # (162.0, 180.0, 198.0)

    results = {
        "capacity_symmetric_sweep": [],
        "capacity_asymmetric_cases": [],
        "demand_forecast_sweep": []
    }

    print("--- 1. Running Capacity Tolerance Symmetric Sweep (Fixed Demand (162, 180, 198)) ---")
    for d in DELTAS:
        b = dict(base_benchmark)
        b["K_tfn"] = fixed_k_tfn

        # 1. Pessimistic (alpha=0, gamma=0) -> Necessity criterion
        res_p = solve_exact_cpp(b, 0.0, 0.0, left_spread=d, right_spread=d)
        T_p = ceil(Fraction(198, 1) / Fraction(str(1 - d)))
        rn_p, rp_p = compute_measures(b, res_p["best_chromosome"], d, d)

        # 2. Neutral (alpha=0.5, gamma=0.5)
        res_n = solve_exact_cpp(b, 0.5, 0.5, left_spread=d, right_spread=d)
        T_n = 180
        rn_n, rp_n = compute_measures(b, res_n["best_chromosome"], d, d)

        # 3. Optimistic (alpha=0, gamma=1) -> Possibility criterion
        res_o = solve_exact_cpp(b, 0.0, 1.0, left_spread=d, right_spread=d)
        T_o = ceil(Fraction(162, 1) / Fraction(str(1 + d)))
        rn_o, rp_o = compute_measures(b, res_o["best_chromosome"], d, d)

        entry = {
            "delta": d,
            "pessimistic": {
                "T": T_p,
                "cost": res_p["best_cost"],
                "fleet": describe_fleet(b, res_p["best_chromosome"]),
                "R_hurwicz": res_p["best_reliability"],
                "R_necessity": rn_p,
                "R_possibility": rp_p,
                "design": res_p["best_chromosome"],
                "nodes": res_p["nodes"],
                "time": res_p["runtime_seconds"],
            },
            "neutral": {
                "T": T_n,
                "cost": res_n["best_cost"],
                "fleet": describe_fleet(b, res_n["best_chromosome"]),
                "R_hurwicz": res_n["best_reliability"],
                "R_necessity": rn_n,
                "R_possibility": rp_n,
                "design": res_n["best_chromosome"],
                "nodes": res_n["nodes"],
                "time": res_n["runtime_seconds"],
            },
            "optimistic": {
                "T": T_o,
                "cost": res_o["best_cost"],
                "fleet": describe_fleet(b, res_o["best_chromosome"]),
                "R_hurwicz": res_o["best_reliability"],
                "R_necessity": rn_o,
                "R_possibility": rp_o,
                "design": res_o["best_chromosome"],
                "nodes": res_o["nodes"],
                "time": res_o["runtime_seconds"],
            }
        }
        results["capacity_symmetric_sweep"].append(entry)
        print(f"delta={d:.2f} | Pess: Cost=${entry['pessimistic']['cost']}M (T={T_p}) | "
              f"Neut: Cost=${entry['neutral']['cost']}M (T={T_n}) | "
              f"Opt: Cost=${entry['optimistic']['cost']}M (T={T_o})")

    print("\n--- 2. Running Capacity Asymmetric Seasonal Cases (Fixed Demand) ---")
    for case in ASYMMETRIC_CASES:
        dl = case["delta_l"]
        dr = case["delta_r"]
        b = dict(base_benchmark)
        b["K_tfn"] = fixed_k_tfn

        res_p = solve_exact_cpp(b, 0.0, 0.0, left_spread=dl, right_spread=dr)
        T_p = ceil(Fraction(198, 1) / Fraction(str(1 - dl)))
        rn_p, rp_p = compute_measures(b, res_p["best_chromosome"], dl, dr)

        res_o = solve_exact_cpp(b, 0.0, 1.0, left_spread=dl, right_spread=dr)
        T_o = ceil(Fraction(162, 1) / Fraction(str(1 + dr)))
        rn_o, rp_o = compute_measures(b, res_o["best_chromosome"], dl, dr)

        # Modal cell (0.5, 0.5): kappa = k2 / m with m = 1 - dl/2 + (dl + dr)/4 (Equation for m)
        res_n = solve_exact_cpp(b, 0.5, 0.5, left_spread=dl, right_spread=dr)
        m_n = 1 - Fraction(str(dl)) / 2 + (Fraction(str(dl)) + Fraction(str(dr))) / 4
        T_n = ceil(Fraction(180, 1) / m_n)
        rn_n, rp_n = compute_measures(b, res_n["best_chromosome"], dl, dr)

        results["capacity_asymmetric_cases"].append({
            "label": case["label"],
            "delta_l": dl,
            "delta_r": dr,
            "pessimistic": {
                "T": T_p,
                "cost": res_p["best_cost"],
                "fleet": describe_fleet(b, res_p["best_chromosome"]),
                "R_hurwicz": res_p["best_reliability"],
                "R_necessity": rn_p,
                "R_possibility": rp_p,
                "design": res_p["best_chromosome"],
            },
            "optimistic": {
                "T": T_o,
                "cost": res_o["best_cost"],
                "fleet": describe_fleet(b, res_o["best_chromosome"]),
                "R_hurwicz": res_o["best_reliability"],
                "R_necessity": rn_o,
                "R_possibility": rp_o,
                "design": res_o["best_chromosome"],
            },
            "neutral": {
                "T": T_n,
                "cost": res_n["best_cost"],
                "fleet": describe_fleet(b, res_n["best_chromosome"]),
                "R_hurwicz": res_n["best_reliability"],
                "R_necessity": rn_n,
                "R_possibility": rp_n,
                "design": res_n["best_chromosome"],
            }
        })
        print(f"{case['label']} | Pess: Cost=${res_p['best_cost']}M (T={T_p}) | "
              f"Neut: Cost=${res_n['best_cost']}M (T={T_n}) | Opt: Cost=${res_o['best_cost']}M (T={T_o})")

    print("\n--- 3. Running Demand Forecast Sweep (Capacity Spread Fixed at 0.10) ---")
    for eta in ETAS:
        b = dict(base_benchmark)
        k1 = round(180.0 * (1.0 - eta), 3)
        k3 = round(180.0 * (1.0 + eta), 3)
        b["K_tfn"] = (k1, 180.0, k3)

        res_p = solve_exact_cpp(b, 0.0, 0.0, left_spread=0.10, right_spread=0.10)
        T_p = ceil(Fraction(str(k3)) / Fraction("0.9"))
        res_n = solve_exact_cpp(b, 0.5, 0.5, left_spread=0.10, right_spread=0.10)
        T_n = 180
        res_o = solve_exact_cpp(b, 0.0, 1.0, left_spread=0.10, right_spread=0.10)
        T_o = ceil(Fraction(str(k1)) / Fraction("1.1"))

        entry = {
            "eta": eta,
            "demand_tfn": [k1, 180.0, k3],
            "pessimistic_cost": res_p["best_cost"],
            "pessimistic_T": T_p,
            "neutral_cost": res_n["best_cost"],
            "neutral_T": T_n,
            "optimistic_cost": res_o["best_cost"],
            "optimistic_T": T_o,
        }
        results["demand_forecast_sweep"].append(entry)

    # Save JSON artifact
    out_json = ROOT / "results" / "generated" / "delta_sensitivity_v7.json"
    out_json.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {out_json}")

    # Plot publication figure
    plot_delta_sensitivity(results)


def plot_delta_sensitivity(results: dict):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)

    # Panel (a): Capacity tolerance sensitivity (symmetric + asymmetric)
    deltas = [item["delta"] for item in results["capacity_symmetric_sweep"]]
    cost_pess = [item["pessimistic"]["cost"] for item in results["capacity_symmetric_sweep"]]
    cost_neut = [item["neutral"]["cost"] for item in results["capacity_symmetric_sweep"]]
    cost_opt = [item["optimistic"]["cost"] for item in results["capacity_symmetric_sweep"]]

    ax1.plot(deltas, cost_pess, "s-", color="#B71C1C", linewidth=2.0, markersize=6,
             label=r"Necessity cell ($\alpha=0, \gamma=0$)")
    ax1.plot(deltas, cost_neut, "o--", color="#1E88E5", linewidth=2.0, markersize=6,
             label=r"Modal cell ($\alpha=0.5, \gamma=0.5$)")
    ax1.plot(deltas, cost_opt, "^-", color="#2E7D32", linewidth=2.0, markersize=6,
             label=r"Optimistic boundary cell ($\alpha=0, \gamma=1$)")

    # Plot asymmetric regime markers dynamically
    # Asymmetric regimes: necessity-cell cost, plotted at x = delta_l
    asym_markers = [
        ("Summer", 0, (8, -14)),
        ("Extreme heat", 1, (-62, 6)),
        ("Winter", 2, (-14, 12)),
        ("Cool climate", 3, (6, -16)),
    ]
    for short, idx, xytext in asym_markers:
        case = results["capacity_asymmetric_cases"][idx]
        dl, cost = case["delta_l"], case["pessimistic"]["cost"]
        ax1.scatter([dl], [cost], color="#6A1B9A", s=55, zorder=5, marker="D")
        ax1.annotate(f"{short} ({case['delta_l']:.2f}, {case['delta_r']:.2f}): ${cost}M", (dl, cost),
                     textcoords="offset points", xytext=xytext, fontsize=7.0, color="#4A148C")

    ax1.axvline(0.10, color="#546E7A", linestyle=":", alpha=0.8, lw=1.2)
    ax1.annotate(r"Baseline $\delta=0.10$", (0.10, 130), textcoords="offset points",
                 xytext=(-65, -30), fontsize=8.0, color="#37474F",
                 arrowprops=dict(arrowstyle="->", color="#546E7A", lw=0.8))

    ax1.set_xlabel(r"Epistemic Capacity Tolerance / Derating $\delta$", fontsize=10.0)
    ax1.set_ylabel(r"Optimal Capital Expenditure $Z^*(\delta)$ (M\$)", fontsize=10.0)
    ax1.set_title(r"(a) Sensitivity to Capacity Derating Half-Width $\delta$", fontsize=10.5, pad=8)
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper left", fontsize=8.2, framealpha=0.95)
    ax1.set_xlim(0.01, 0.21)
    ax1.set_ylim(110, 290)

    # Panel (b): Demand Forecast Uncertainty Sensitivity (eta in [0.02, 0.20])
    etas = [item["eta"] for item in results["demand_forecast_sweep"]]
    d_pess = [item["pessimistic_cost"] for item in results["demand_forecast_sweep"]]
    d_neut = [item["neutral_cost"] for item in results["demand_forecast_sweep"]]
    d_opt = [item["optimistic_cost"] for item in results["demand_forecast_sweep"]]

    ax2.plot(etas, d_pess, "s-", color="#B71C1C", linewidth=2.0, markersize=6,
             label=r"Necessity cell ($\alpha=0, \gamma=0$)")
    ax2.plot(etas, d_neut, "o--", color="#1E88E5", linewidth=2.0, markersize=6,
             label=r"Modal cell ($\alpha=0.5, \gamma=0.5$)")
    ax2.plot(etas, d_opt, "^-", color="#2E7D32", linewidth=2.0, markersize=6,
             label=r"Optimistic boundary cell ($\alpha=0, \gamma=1$)")

    ax2.axvline(0.10, color="#546E7A", linestyle=":", alpha=0.8, lw=1.2)
    ax2.annotate(r"Baseline $\eta=0.10$", (0.10, 130), textcoords="offset points",
                 xytext=(-65, -30), fontsize=8.0, color="#37474F",
                 arrowprops=dict(arrowstyle="->", color="#546E7A", lw=0.8))

    ax2.set_xlabel(r"Peak Demand Forecast Uncertainty Half-Width $\eta$", fontsize=10.0)
    ax2.set_ylabel(r"Optimal Capital Expenditure $Z^*(\eta)$ (M\$)", fontsize=10.0)
    ax2.set_title(r"(b) Sensitivity to Peak Demand Uncertainty $\eta$", fontsize=10.5, pad=8)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper left", fontsize=8.2, framealpha=0.95)
    ax2.set_xlim(0.01, 0.21)
    ax2.set_ylim(110, 290)

    fig.tight_layout()
    fig_pdf = ROOT / "figures" / "fig_delta_sensitivity.pdf"
    fig_png = ROOT / "figures" / "fig_delta_sensitivity.png"
    fig.savefig(fig_pdf, bbox_inches="tight")
    fig.savefig(fig_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fig_pdf} and {fig_png}")


if __name__ == "__main__":
    run_delta_sensitivity()
