"""Validate every native Pareto-DP result with the Python reference kernel."""

from __future__ import annotations

import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.solver import evaluate_design


def main() -> None:
    source = ROOT / "results" / "generated" / "fyffe_pareto_dp.csv"
    destination = ROOT / "results" / "generated" / "fyffe_exact_audit.json"
    benchmark = BENCHMARKS["fyffe_adapted"]
    rows = list(csv.DictReader(source.open(encoding="utf-8")))
    if len(rows) != 30:
        raise AssertionError(f"Expected 30 cells, found {len(rows)}")

    seen: set[tuple[float, float]] = set()
    audit: list[dict] = []
    for row in rows:
        alpha, gamma = float(row["alpha"]), float(row["gamma"])
        key = (alpha, gamma)
        if key in seen:
            raise AssertionError(f"Duplicate cell {key}")
        seen.add(key)
        chromosome = [int(value) for value in row["chromosome"].split("-")]
        if len(chromosome) != benchmark["N"]:
            raise AssertionError(f"Invalid chromosome length at {key}")
        checked = evaluate_design(chromosome, benchmark, alpha, gamma)
        native_cost = int(row["best_cost"])
        native_reliability = float(row["best_reliability"])
        if checked["cost"] != native_cost:
            raise AssertionError(f"Cost mismatch at {key}")
        if abs(checked["reliability"] - native_reliability) > 1e-10:
            raise AssertionError(f"Reliability mismatch at {key}")
        if not checked["feasible"]:
            raise AssertionError(f"Infeasible native result at {key}")
        audit.append({
            "alpha": alpha,
            "gamma": gamma,
            "best_cost": native_cost,
            "best_reliability_native": native_reliability,
            "best_reliability_python": checked["reliability"],
            "absolute_difference": abs(checked["reliability"] - native_reliability),
            "chromosome": chromosome,
            "threshold_scaled_python": checked["threshold_scaled"],
            "scale_python": checked["scale"],
            "peak_native_states": int(row["peak_states"]),
            "runtime_seconds_native": float(row["seconds"]),
            "passed": True,
        })

    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "benchmark": "fyffe_adapted",
        "cells_expected": 30,
        "cells_verified": len(audit),
        "all_passed": all(item["passed"] for item in audit),
        "objective": ["min_cost", "max_reliability_at_minimum_cost"],
        "exactness_basis": (
            "Cost-bounded enumeration of every DP state up to a known feasible cost; "
            "within each layer and exact cost only first-order stochastically dominated "
            "capped-capacity PMFs are removed. A feasible cap excludes all higher costs "
            "from improving the primary objective."
        ),
        "independent_check": "Every reported chromosome was recomputed by the Python rational-scaling capped-PMF kernel.",
        "results": audit,
    }
    destination.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Verified {len(audit)} Fyffe cells; wrote {destination}")


if __name__ == "__main__":
    main()
