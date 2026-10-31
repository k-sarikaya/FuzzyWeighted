"""Full Pharma oracle audit; intentionally separate from the fast unit suite."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.oracle import vectorized_exhaustive_optimum
from fuzzy_wkon.solver import solve_exact


def main() -> int:
    benchmark = BENCHMARKS["pharma_adapted"]
    alphas = (0.0, 0.3, 0.5, 0.7, 0.9, 1.0)
    gammas = (0.0, 0.3, 0.5, 0.7, 1.0)
    audit = []
    for alpha in alphas:
        for gamma in gammas:
            exact = solve_exact(benchmark, alpha, gamma)
            oracle = vectorized_exhaustive_optimum(benchmark, alpha, gamma)
            match = (
                exact["feasible"] == oracle["feasible"]
                and exact["best_cost"] == oracle["cost"]
                and (not exact["feasible"] or abs(exact["best_reliability"] - oracle["reliability"]) <= 1e-12)
            )
            row = {"alpha": alpha, "gamma": gamma, "match": match, "solver": exact, "oracle": oracle}
            audit.append(row)
            print(json.dumps({"alpha": alpha, "gamma": gamma, "match": match, "cost": exact["best_cost"]}))
            if not match:
                return 1
    out = ROOT / "results" / "generated" / "pharma_oracle_audit.json"
    out.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(f"Verified {len(audit)} Pharma cells; wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
