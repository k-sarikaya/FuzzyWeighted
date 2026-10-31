"""Enumerate all 3^10 = 59,049 RTS-79 designs to compute:
1. True Cost-Necessity (Z vs R_N) Pareto frontier (min Z s.t. R_N >= lambda, max R_N s.t. Z <= B).
2. Compare with Hurwicz menu designs (D1..D7).
3. Find the exact counterexample design ($171M) noted by Reviewer R2.
4. R_target sensitivity (0.90, 0.95, 0.98) showing where redundancy becomes active.
"""

import sys
import json
import itertools
from collections import Counter
from fractions import Fraction
from math import lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.fuzzy import TFN
from fuzzy_wkon.membership import membership_profile, saturation_level
from fuzzy_wkon.reliability import measure_reliability, static_reliability

SPREAD = 0.10
MODAL_DEMAND = 180          # crisp adequacy is judged at the modal demand k2
R_TARGET = 0.95
EPS = 1e-12
B = BENCHMARKS["power_rts_adapted"]
K_TFN = TFN(*B["K_tfn"])
scale = 1

def modal_options(benchmark: dict, chromosome: list[int]) -> tuple[list[tuple[float, int]], int]:
    selected: list[tuple[float, int]] = []
    for choices, choice_index in zip(benchmark["choices"], chromosome):
        reliability, _cost, weight = choices[choice_index]
        if choice_index == 0:
            continue
        selected.append((float(reliability), int(weight)))
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
    desc = " + ".join(parts) if parts else "Empty fleet"
    return desc, total_mw, total_cost

cap = max(saturation_level(K_TFN, "pos", SPREAD, scale), saturation_level(K_TFN, "nec", SPREAD, scale))
nec_profile = membership_profile(K_TFN, "nec", SPREAD, SPREAD, scale, cap)
pos_profile = membership_profile(K_TFN, "pos", SPREAD, SPREAD, scale, cap)

print(f"Capacity cap: {cap}")
print(f"Necessity profile len: {len(nec_profile)}, Possibility profile len: {len(pos_profile)}")

# Enumerate all 3^10 combinations
choices_list = B["choices"]
num_bays = len(choices_list)
all_combos = list(itertools.product(range(3), repeat=num_bays))
print(f"Total combinatorial designs: {len(all_combos)}")

catalog_results = []
for idx, combo in enumerate(all_combos):
    cost = sum(choices_list[i][combo[i]][1] for i in range(num_bays))
    total_w = sum(choices_list[i][combo[i]][2] for i in range(num_bays))
    options, _ = modal_options(B, list(combo))
    
    r_n = measure_reliability(options, nec_profile) if options else 0.0
    r_p = measure_reliability(options, pos_profile) if options else 0.0
    r_modal = static_reliability(options, MODAL_DEMAND) if options else 0.0

    catalog_results.append({
        "combo": list(combo),
        "cost": cost,
        "total_w": total_w,
        "r_n": round(r_n, 6),
        "r_p": round(r_p, 6),
        "r_modal": round(r_modal, 6),
        "adequate": r_modal >= R_TARGET - EPS,
    })

print("Enumeration finished!")

# Find designs with cost 171
cost_171 = [d for d in catalog_results if d["cost"] == 171]
print(f"\nFound {len(cost_171)} designs with cost = 171:")
for d in cost_171:
    desc, mw, c = describe_fleet(B, d["combo"])
    print(f"   Fleet: {desc:35s} | MW: {mw:3d} | Cost: {c:3d} | R_N: {d['r_n']:.4f} | R_Pi: {d['r_p']:.4f}")

# Compute True Pareto Frontier (Cost vs R_N)
# Non-dominated designs: for increasing cost, r_n must be strictly greater than max seen so far
pareto_frontier = []
catalog_sorted = sorted(catalog_results, key=lambda x: (x["cost"], -x["r_n"]))

max_rn_so_far = -1.0
for d in catalog_sorted:
    if d["r_n"] > max_rn_so_far + 1e-6:
        desc, mw, c = describe_fleet(B, d["combo"])
        d_copy = dict(d)
        d_copy["fleet"] = desc
        pareto_frontier.append(d_copy)
        max_rn_so_far = d["r_n"]

print(f"\nTrue Pareto Frontier size (min Cost, max R_N): {len(pareto_frontier)}")
for p in pareto_frontier:
    print(f"Cost: {p['cost']:3d} | R_N: {p['r_n']:.4f} | R_Pi: {p['r_p']:.4f} | MW: {p['total_w']:3d} | Fleet: {p['fleet']}")

# Compare with Hurwicz menu designs D1..D7
menu_path = ROOT / "results" / "generated" / "policy_menu_v8.json"
menu_data = None
if menu_path.exists():
    menu_data = json.loads(menu_path.read_text(encoding="utf-8"))
    print("\nHurwicz Policy Menu (D1..D7):")
    for d in menu_data["menu"]:
        print(f"   {d['label']}: Cost: {d['cost']:3d} | R_N: {d['R_N']:.4f} | R_Pi: {d['R_Pi']:.4f} | Fleet: {d['fleet']}")

# Adequacy-constrained frontier: min Z s.t. R_N >= lambda and R_modal >= R_target
adequate_frontier = []
max_rn_so_far = -1.0
for d in catalog_sorted:
    if d["adequate"] and d["r_n"] > max_rn_so_far + 1e-6:
        desc, mw, c = describe_fleet(B, d["combo"])
        d_copy = dict(d)
        d_copy["fleet"] = desc
        adequate_frontier.append(d_copy)
        max_rn_so_far = d["r_n"]
cheapest_adequate = min((d for d in catalog_results if d["adequate"]), key=lambda d: d["cost"])

by_combo = {tuple(d["combo"]): d for d in catalog_results}
menu_modal = []
if menu_data is not None:
    for m in menu_data["menu"]:
        d = by_combo[tuple(m["chromosome"])]
        menu_modal.append({"label": m["label"], "cost": m["cost"], "total_w": d["total_w"],
                           "r_modal": d["r_modal"], "adequate": d["adequate"]})
menu_combos = {tuple(m["chromosome"]) for m in menu_data["menu"]} if menu_data else set()

print("\nAdequacy-constrained frontier (R_modal >= 0.95):")
for p in adequate_frontier:
    print(f"Cost: {p['cost']:3d} | R_N: {p['r_n']:.4f} | R_modal: {p['r_modal']:.4f} | Fleet: {p['fleet']}")
print("Cheapest adequate design cost:", cheapest_adequate["cost"])
print("Menu R_modal:", [(m["label"], m["r_modal"]) for m in menu_modal])

# LaTeX table of the adequacy-constrained frontier between D6 and D7
rows = []
for p in adequate_frontier:
    if p["cost"] > 203:
        break
    label = next((m["label"] for m in menu_data["menu"] if tuple(m["chromosome"]) == tuple(p["combo"])), "--") if menu_data else "--"
    fleet = p["fleet"].replace("×", r" $\times$ ")
    rows.append(f"{p['cost']} & {fleet} & {p['total_w']} & {p['r_modal']:.4f} & {p['r_n']:.4f} & {p['r_p']:.4f} & {label} \\\\")
table = (
    "\\begin{table}[!htbp]\n\\centering\n"
    "\\caption{Adequacy-constrained cost--necessity frontier of the IEEE RTS-79 catalog: designs that minimize cost "
    "subject to $R_N(\\mathbf{x}) \\ge r$ and crisp adequacy $R_{\\mathrm{modal}}(\\mathbf{x}) \\ge 0.95$ "
    "at the modal demand of $180$~MW, for every $r$ up to the necessity design D7 (exhaustive enumeration of all "
    "$59{,}049$ designs). The last column names the menu design (Table~\\ref{tab:policy_menu}) when the "
    "frontier design belongs to it.}\n"
    "\\label{tab:adequate_frontier}\n\\footnotesize\n\\begin{tabular}{rp{4.6cm}rcccc}\n\\hline\n"
    "Cost (M\\$) & Fleet & MW & $R_{\\mathrm{modal}}$ & $R_N$ & $R_\\Pi$ & Menu\\\\\n\\hline\n"
    + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n\\end{table}\n")
(ROOT / "manuscript" / "table_adequate_frontier.tex").write_text(table, encoding="utf-8")

out_path = ROOT / "results" / "generated" / "rts_true_frontier_v9.json"
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w", encoding="utf-8") as f:
    json.dump({
        "total_evaluated": len(catalog_results),
        "pareto_frontier": pareto_frontier,
        "cost_171_designs": cost_171,
        "adequacy": {"modal_demand": MODAL_DEMAND, "R_target": R_TARGET},
        "adequate_frontier": adequate_frontier,
        "cheapest_adequate": cheapest_adequate,
        "menu_modal": menu_modal,
        "unconstrained_frontier_adequate_costs": [p["cost"] for p in pareto_frontier
                                                  if by_combo[tuple(p["combo"])]["adequate"]],
        "adequate_frontier_not_in_menu": [p["cost"] for p in adequate_frontier
                                          if tuple(p["combo"]) not in menu_combos and p["cost"] <= 203],
    }, f, indent=2)

print(f"\nSaved results to {out_path}")
