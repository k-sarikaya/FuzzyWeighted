"""Complete 30-cell enumeration of IEEE RTS-79:
1. Baseline common spread (0.10, 0.10)
2. Unit-specific spreads (Table 6 RTS-79)
3. Matched common spread (0.06, 0.02) - isolates heterogeneity from uncertainty level

Resolves Reviewer concerns R2-M4 and R3-M4.
"""

import sys
import json
import itertools
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS

B = BENCHMARKS["power_rts_adapted"]
ALPHAS = (0.0, 0.3, 0.5, 0.7, 0.9, 1.0)
GAMMAS = (0.0, 0.3, 0.5, 0.7, 1.0)

COMMON_10 = {mw: (0.10, 0.10) for mw in (12, 20, 50, 76, 100, 155)}
UNIT_SPEC = {12: (0.05, 0.02), 20: (0.15, 0.05), 50: (0.10, 0.0), 76: (0.05, 0.02), 100: (0.05, 0.02), 155: (0.05, 0.02)}
MATCHED_AVG = {mw: (0.06, 0.02) for mw in (12, 20, 50, 76, 100, 155)}

# Precompute designs
all_designs = list(itertools.product(*[range(len(c)) for c in B["choices"]]))
design_data = []
for d in all_designs:
    units = [B["choices"][i][c] for i, c in enumerate(d) if c]
    cost = sum(u[1] for u in units)
    mws = [u[2] for u in units]
    probs = [F(str(u[0])) for u in units]
    design_data.append((cost, mws, probs, d, units))

# Sort designs by cost ascending
design_data.sort(key=lambda x: x[0])

def solve_cell(alpha: float, gamma: float, spreads: dict):
    a, g = F(str(alpha)), F(str(gamma))
    k1, k2, k3 = (F(str(x)) for x in B["K_tfn"])
    lower, upper = k1 + a * (k2 - k1), k3 - a * (k3 - k2)
    threshold = upper - g * (upper - lower)

    def cap_func(mw: int) -> F:
        left, right = (F(str(x)) for x in spreads[mw])
        w = F(mw)
        lo = w * (1 - left) + a * left * w
        hi = w * (1 + right) - a * right * w
        return lo + g * (hi - lo)

    # Memoize capacities for the 6 MW ratings
    memo_cap = {mw: cap_func(mw) for mw in (12, 20, 50, 76, 100, 155)}
    target = F(str(B["R_target"]))

    best = None
    for cost, mws, probs, d, units in design_data:
        if best is not None and cost > best[0]:
            break  # because design_data is sorted by cost
        caps = [memo_cap[mw] for mw in mws]
        if sum(caps) < threshold:
            continue
        
        # Exact reliability
        rel = F(0)
        n_u = len(units)
        for states in itertools.product((0, 1), repeat=n_u):
            if sum(c for c, s in zip(caps, states) if s) >= threshold:
                p = F(1)
                for pu, s in zip(probs, states):
                    p *= pu if s else 1 - pu
                rel += p
        
        if rel < target:
            continue
        
        if best is None or cost < best[0] or (cost == best[0] and rel > best[1]):
            best = (cost, rel, d)
            
    return best

print("Running 30-cell comparison: Common-10% vs Unit-Specific vs Matched-Common (0.06, 0.02)...")
results = []
diff_count_unit_vs_common = 0
diff_count_unit_vs_matched = 0

for g in GAMMAS:
    for a in ALPHAS:
        sol_c10 = solve_cell(a, g, COMMON_10)
        sol_unit = solve_cell(a, g, UNIT_SPEC)
        sol_match = solve_cell(a, g, MATCHED_AVG)
        
        cost_c10 = sol_c10[0] if sol_c10 else None
        cost_unit = sol_unit[0] if sol_unit else None
        cost_match = sol_match[0] if sol_match else None
        
        d_c10 = sol_c10[2] if sol_c10 else None
        d_unit = sol_unit[2] if sol_unit else None
        d_match = sol_match[2] if sol_match else None
        
        c_vs_u = (d_c10 != d_unit)
        u_vs_m = (d_unit != d_match)
        if c_vs_u: diff_count_unit_vs_common += 1
        if u_vs_m: diff_count_unit_vs_matched += 1
        
        row = {
            "alpha": a, "gamma": g,
            "cost_common10": cost_c10, "design_common10": list(d_c10) if d_c10 else None,
            "cost_unit": cost_unit, "design_unit": list(d_unit) if d_unit else None,
            "cost_matched": cost_match, "design_matched": list(d_match) if d_match else None,
            "diff_unit_vs_common": c_vs_u,
            "diff_unit_vs_matched": u_vs_m
        }
        results.append(row)
        print(f"a={a:.1f}, g={g:.1f} | C10: ${cost_c10}M | Unit: ${cost_unit}M | Match: ${cost_match}M | diff(U,C10): {c_vs_u} | diff(U,Match): {u_vs_m}")

print(f"\nTotal cells where Unit-Specific differs from Common-10%: {diff_count_unit_vs_common} / 30")
print(f"Total cells where Unit-Specific differs from Matched-Common: {diff_count_unit_vs_matched} / 30")

out_file = ROOT / "results" / "generated" / "rts_unit_specific_matched_v9.json"
out_file.parent.mkdir(parents=True, exist_ok=True)
out_file.write_text(json.dumps({
    "diff_unit_vs_common": diff_count_unit_vs_common,
    "diff_unit_vs_matched": diff_count_unit_vs_matched,
    "cells": results
}, indent=2), encoding="utf-8")
print(f"Saved results to {out_file}")
