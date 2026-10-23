"""Sensitivity to target reliability R_target in IEEE RTS-79 fleet.

Resolves Reviewer concern R2-M1(b) and R2-M1(c):
- Tests R_target in {0.90, 0.92, 0.95, 0.97, 0.98, 0.99}.
- Evaluates whether redundancy becomes active (i.e. units in excess of minimum series capacity).
- Identifies how many single-unit losses each design can withstand while still meeting T.
"""

import sys
import json
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.fuzzy import TFN
from fuzzy_wkon.cpp_backend import solve_at_integer_threshold_cpp
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "campaign", ROOT / "experiments" / "run_power_rts_campaign_v7.py")
campaign = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(campaign)

SPREAD = 0.10
B_ORIG = BENCHMARKS["power_rts_adapted"]
K_TFN = B_ORIG["K_tfn"]

# Test for modal condition (T=180), pessimistic (T=220, cell 0,0), and optimistic (T=148, cell 0,1)
T_TESTS = [
    {"label": "Optimistic boundary (0,1)", "T": 148},
    {"label": "Modal (0.5,0.5)", "T": 180},
    {"label": "Necessity (0,0)", "T": 220}
]

R_TARGETS = [round(0.90 + 0.01 * i, 2) for i in range(10)]   # 0.90, 0.91, ..., 0.99

results = []

for case in T_TESTS:
    T = case["T"]
    case_label = case["label"]
    print(f"\n--- Testing {case_label} with T = {T} MW ---")
    for r_target in R_TARGETS:
        b = deepcopy(B_ORIG)
        b["R_target"] = r_target
        out = solve_at_integer_threshold_cpp(b, T, time_limit=60.0)
        feasible = out.get("status") == "completed" and out.get("feasible")
        
        if feasible:
            chromo = out["best_chromosome"]
            desc, mw, cost = campaign.describe_fleet(b, chromo)
            r_hurwicz = out["best_reliability"]
            
            # Check unit loss tolerance: can the fleet lose any unit and still have MW >= T?
            # List installed units
            installed_units = []
            for bay, c in enumerate(chromo):
                if c != 0:
                    installed_units.append(b["choices"][bay][c][2]) # MW
            
            # Check minimum capacity after losing 1 unit
            min_after_1_loss = mw - max(installed_units) if installed_units else 0
            tolerates_single_loss = (min_after_1_loss >= T)
            
            # Count redundancy: total installed MW minus T
            reserve_margin = mw - T
            
            row = {
                "case": case_label,
                "T": T,
                "R_target": r_target,
                "feasible": True,
                "cost": cost,
                "total_mw": mw,
                "reserve_mw": reserve_margin,
                "fleet": desc,
                "R_sys": round(r_hurwicz, 5),
                "tolerates_single_loss": tolerates_single_loss,
                "min_mw_after_1_loss": min_after_1_loss,
                "nodes": out.get("nodes"),
                "time": round(out.get("runtime_seconds", 0), 3)
            }
        else:
            row = {
                "case": case_label,
                "T": T,
                "R_target": r_target,
                "feasible": False,
                "cost": None,
                "fleet": "Infeasible"
            }
        results.append(row)
        if row["feasible"]:
            print(f"R_target={r_target:.2f} | Cost=${row['cost']:3d}M | MW={row['total_mw']:3d} (T={T}) | R_sys={row['R_sys']:.4f} | 1-Loss Tolerant: {row['tolerates_single_loss']} | Fleet: {row['fleet']}")
        else:
            print(f"R_target={r_target:.2f} | INFEASIBLE")

out_file = ROOT / "results" / "generated" / "rts_rtarget_sensitivity_v9.json"
out_file.parent.mkdir(parents=True, exist_ok=True)
out_file.write_text(json.dumps(results, indent=2), encoding="utf-8")
print(f"\nSaved R_target sensitivity results to {out_file}")


def fleet_tex(desc: str) -> str:
    return desc.replace("×", r" $\times$ ")


def interval_rows(rows: list[dict]) -> list[tuple[float, float, dict]]:
    """Group consecutive R_target values with the same optimal design."""
    groups = []
    for r in rows:
        if groups and groups[-1][2]["fleet"] == r["fleet"] and groups[-1][2]["cost"] == r["cost"]:
            groups[-1] = (groups[-1][0], r["R_target"], groups[-1][2])
        else:
            groups.append((r["R_target"], r["R_target"], r))
    return groups


lines = []
for case in T_TESTS:
    rows = [r for r in results if r["case"] == case["label"] and r["feasible"]]
    groups = interval_rows(rows)
    for j, (lo, hi, r) in enumerate(groups):
        regime = (f"\\multirow{{{len(groups)}}}{{*}}{{\\shortstack[l]{{{case['label'].split(' (')[0]}\\\\"
                  f"$({case['label'].split('(')[1][:-1]})$, $T={case['T']}$~MW}}}}") if j == 0 else ""
        rng = f"{lo:.2f}" if lo == hi else f"{lo:.2f}--{hi:.2f}"
        n1 = (f"Yes ({r['min_mw_after_1_loss']} MW)" if r["tolerates_single_loss"]
              else f"No ({r['min_mw_after_1_loss']} MW)")
        lines.append(f"{regime} & {rng} & {fleet_tex(r['fleet'])} & {r['total_mw']} & {r['cost']} & "
                     f"{r['R_sys']:.4f} & {n1}\\\\")
    lines.append("\\hline")
table = (
    "\\begin{table}[!htbp]\n\\centering\n"
    "\\caption{Sensitivity of the optimal fleet to the reliability target $R_{\\mathrm{target}} \\in \\{0.90, 0.91, "
    "\\dots, 0.99\\}$ on the IEEE RTS-79 fleet at three effective thresholds. Consecutive targets with the same "
    "optimal design are merged. $R_{\\mathrm{sys}}$ is the reliability at the effective threshold $T$. The last column states "
    "whether the working capacity still reaches $T$ after the loss of the largest installed unit, and gives that "
    "remaining capacity.}\n"
    "\\label{tab:rts_rtarget_sensitivity}\n\\scriptsize\n\\setlength{\\tabcolsep}{1.5pt}\n"
    "\\begin{tabular}{lcp{4.4cm}cccc}\n\\hline\n"
    "Regime & $R_{\\mathrm{target}}$ & Optimal fleet & MW & Cost (M\\$) & $R_{\\mathrm{sys}}$ & Largest-unit loss\\\\\n\\hline\n"
    + "\n".join(lines) + "\n\\end{tabular}\n\\end{table}\n")
(ROOT / "manuscript" / "table_rts_rtarget_sensitivity.tex").write_text(table, encoding="utf-8")
print("Wrote manuscript/table_rts_rtarget_sensitivity.tex")
