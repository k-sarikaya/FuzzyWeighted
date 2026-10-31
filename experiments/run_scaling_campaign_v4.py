"""Bounded scaling campaign, C++ backend plus the Proposition-E reduction.

Two changes relative to the V2 campaign, both of which only affect runtime:

* the search runs in the validated C++ port rather than the Python kernel;
* the instance is reduced to its modal integer form through the effective
  threshold ``T = ceil(kappa)``, which is exact because the modal capacity sum
  is an integer.

The reduction is *not* assumed: `--verify-reduction` re-solves each completed
instance in the unreduced formulation and requires identical cost, reliability
and chromosome before the row is accepted.

Results are appended to disk after every instance, so an interrupted run keeps
its completed observations. Timeouts are retained as censored observations and
are never reported as exact optima.

    PYTHONPATH=src python3 experiments/run_scaling_campaign_v4.py \
        --timeout 600 --sizes 8 10 12 14 16 18 20
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from fractions import Fraction
from math import ceil
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.cpp_backend import (  # noqa: E402
    ensure_binary,
    solve_at_integer_threshold_cpp,
    solve_exact_cpp,
)
from fuzzy_wkon.scaling import PROFILES, synthetic_benchmark  # noqa: E402

ALPHA, GAMMA, SPREAD = 0.3, 0.7, 0.1
DESTINATION = ROOT / "results" / "generated" / "scaling_campaign_v4.json"


def effective_threshold(threshold_tfn, alpha: float, gamma: float, spread: float) -> int:
    """T = ceil(kappa); see docs/V4_THEORY_NOTES.md, Proposition E."""

    a, g, d = Fraction(str(alpha)), Fraction(str(gamma)), Fraction(str(spread))
    k1, k2, k3 = (Fraction(str(value)) for value in threshold_tfn)
    cut_lower = k1 + a * (k2 - k1)
    cut_upper = k3 - a * (k3 - k2)
    demand = cut_upper - g * (cut_upper - cut_lower)
    return ceil(demand / (1 + d * (1 - a) * (2 * g - 1)))


def solve_one(n: int, profile: str, seed: int, timeout: float, verify: bool) -> dict:
    benchmark = synthetic_benchmark(n, profile, seed)
    threshold = effective_threshold(benchmark["K_tfn"], ALPHA, GAMMA, SPREAD)

    started = perf_counter()
    result = solve_at_integer_threshold_cpp(benchmark, threshold, time_limit=timeout)
    wall = perf_counter() - started

    row = {
        "n": n,
        "profile": profile,
        "seed": seed,
        "alpha": ALPHA,
        "gamma": GAMMA,
        "effective_threshold": threshold,
        "status": result["status"],
        "wall_seconds": wall,
        "nodes": result.get("nodes"),
        "search_space": 3 ** n,
        "backend": "cpp+reduction",
        "generator": benchmark["generator"],
    }
    if result["status"] == "completed":
        row.update({
            "feasible": result.get("feasible"),
            "best_cost": result.get("best_cost"),
            "best_reliability": result.get("best_reliability"),
            "best_chromosome": result.get("best_chromosome"),
            "leaves": result.get("leaves"),
        })
        if verify and result.get("feasible"):
            check_started = perf_counter()
            direct = solve_exact_cpp(benchmark, ALPHA, GAMMA, time_limit=timeout)
            row["reduction_check"] = {
                "status": direct["status"],
                "wall_seconds": perf_counter() - check_started,
                "unreduced_threshold": direct.get("threshold_scaled"),
            }
            if direct["status"] == "completed":
                agrees = (direct.get("best_cost") == result.get("best_cost")
                          and abs(direct.get("best_reliability", 0.0)
                                  - result.get("best_reliability", 0.0)) < 1e-9
                          and direct.get("best_chromosome") == result.get("best_chromosome"))
                row["reduction_check"]["agrees"] = bool(agrees)
                if not agrees:
                    row["reduction_check"]["unreduced"] = {
                        "best_cost": direct.get("best_cost"),
                        "best_reliability": direct.get("best_reliability"),
                        "best_chromosome": direct.get("best_chromosome"),
                    }
    else:
        row["timeout_seconds"] = timeout
    return row


def write(rows: list[dict], timeout: float, sizes: list[int]) -> None:
    DESTINATION.write_text(json.dumps({
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "bounded computational profile, not a general scalability claim",
        "solver": "two-phase exact branch-and-bound, C++ port of the audited Python search",
        "reduction": "modal integer form via T = ceil(kappa) (Proposition E); exact because "
                     "the modal capacity sum is integral",
        "sizes": sizes,
        "profiles": PROFILES,
        "timeout_seconds_per_instance": timeout,
        "runtime_environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "compiler": "g++ -O3 -std=c++17",
        },
        "rows": rows,
    }, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--sizes", nargs="*", type=int,
                        default=[8, 10, 12, 14, 16, 18, 20])
    parser.add_argument("--verify-reduction", action="store_true")
    args = parser.parse_args()

    ensure_binary()
    rows: list[dict] = []
    for n in args.sizes:
        for profile in PROFILES:
            row = solve_one(n, profile, args.seed, args.timeout, args.verify_reduction)
            rows.append(row)
            write(rows, args.timeout, args.sizes)  # checkpoint after every instance
            marker = row["status"]
            if row.get("reduction_check") is not None:
                marker += f" (reduction agrees={row['reduction_check'].get('agrees')})"
            print(f"n={n:>3} {profile:<10} T={row['effective_threshold']:>4} "
                  f"{marker:<32} {row['wall_seconds']:8.2f}s "
                  f"cost={row.get('best_cost')} nodes={row.get('nodes')}", flush=True)
    print(f"\nWrote {len(rows)} observations to {DESTINATION}")


if __name__ == "__main__":
    main()
