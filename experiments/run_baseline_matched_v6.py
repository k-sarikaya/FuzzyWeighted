"""Budget-matched GA comparison with equal-cost polish (round-2 reviewer response).

Round 1 compared a Python GA that stopped on its own stagnation rule after
~1e5 evaluations against a C++ exact solver that ran for 371 s. Two reviewers
ruled that comparison uninformative, and they were right: it could not
distinguish "the certificate buys something" from "the heuristic was starved".

This experiment removes every part of that objection that can be removed
without rewriting the GA in C++, and reports the outcome whatever it is.

Protocol A -- wall-clock matched. Every instance, 5 seeds. The GA receives the
same wall-clock the exact solver used on that instance (floored at 5 s so that
sub-second instances are not degenerate). The GA restarts on stagnation rather
than exiting, so it spends the whole budget.

Protocol B -- evaluation-rich. The two design-profile instances where the
secondary objective was missed or contested, 5 seeds, 1,000,000 evaluations
each. Wall-clock matching favours the compiled solver because the Python GA is
slower per candidate; an evaluation budget is the currency that does not.

Both protocols apply an equal-cost local search to the GA's incumbent -- the
fix a practitioner would apply to a cost-driven search whose secondary
objective is reliability.

    PYTHONPATH=src python3 experiments/run_baseline_matched_v6.py
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from fractions import Fraction
from math import ceil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS  # noqa: E402
from fuzzy_wkon.cpp_backend import solve_at_integer_threshold_cpp  # noqa: E402
from fuzzy_wkon.heuristics import run_ga_matched  # noqa: E402
from fuzzy_wkon.scaling import PROFILES, synthetic_benchmark  # noqa: E402

RESULTS = ROOT / "results" / "generated"
DESTINATION = RESULTS / "baseline_matched_v6.json"
SEEDS = (0, 1, 2, 3, 4)
ALPHA, GAMMA, SPREAD = 0.3, 0.7, 0.1
WALL_FLOOR = 5.0
EVALUATION_BUDGET = 1_000_000
PROTOCOL_B = {("synthetic", 20, "design"), ("synthetic", 14, "design")}


def effective_threshold(threshold_tfn, alpha, gamma, spread) -> int:
    a, g, d = Fraction(str(alpha)), Fraction(str(gamma)), Fraction(str(spread))
    k1, k2, k3 = (Fraction(str(v)) for v in threshold_tfn)
    lower = k1 + a * (k2 - k1)
    upper = k3 - a * (k3 - k2)
    return ceil((upper - g * (upper - lower)) / (1 + d * (1 - a) * (2 * g - 1)))


def crisp(benchmark: dict, threshold: int) -> dict:
    instance = dict(benchmark)
    instance["K_tfn"] = (float(threshold), float(threshold), float(threshold))
    return instance


def summarise(runs: list[dict], exact: dict) -> dict:
    feasible = [r for r in runs if r.get("feasible")]
    cost_hits = [r for r in feasible if r["cost"] == exact["best_cost"]]
    secondary = [r for r in cost_hits
                 if abs(r["reliability"] - exact["best_reliability"]) < 1e-9]
    before = [r for r in cost_hits
              if abs(r.get("reliability_before_polish", r["reliability"])
                     - exact["best_reliability"]) < 1e-9]
    return {
        "cost_hits": len(cost_hits),
        "secondary_hits_after_polish": len(secondary),
        "secondary_hits_before_polish": len(before),
        "seeds": len(SEEDS),
        "worst_reliability_shortfall_at_optimal_cost": (
            max(exact["best_reliability"] - r["reliability"] for r in cost_hits)
            if cost_hits else None),
        "polish_gain_max": max((r.get("polish_gain", 0.0) for r in feasible), default=0.0),
        "mean_evaluations": (sum(r["evaluations"] for r in runs) / len(runs)) if runs else 0,
        "mean_wall_seconds": (sum(r["wall_seconds"] for r in runs) / len(runs)) if runs else 0,
    }


def instances():
    for key, benchmark in BENCHMARKS.items():
        yield ("catalog", key, benchmark)
    for n in (12, 14, 16, 18, 20):
        for profile in PROFILES:
            yield ("synthetic", (n, profile), synthetic_benchmark(n, profile))


def main() -> None:
    rows = []
    for kind, tag, benchmark in instances():
        threshold = effective_threshold(benchmark["K_tfn"], ALPHA, GAMMA, SPREAD)
        exact = solve_at_integer_threshold_cpp(benchmark, threshold, time_limit=1200.0)
        if exact.get("status") != "completed" or not exact.get("feasible"):
            continue

        label = tag if kind == "catalog" else f"n={tag[0]}/{tag[1]}"
        wall = max(exact["runtime_seconds"], WALL_FLOOR)

        row = {
            "label": label,
            "n": benchmark["N"],
            "T": threshold,
            "exact_cost": exact["best_cost"],
            "exact_reliability": exact["best_reliability"],
            "exact_runtime_seconds": exact["runtime_seconds"],
            "protocol_a": {"wall_budget_seconds": wall, "runs": []},
        }

        for seed in SEEDS:
            row["protocol_a"]["runs"].append(
                run_ga_matched(crisp(benchmark, threshold), 0.0, 0.0,
                               left_spread=0.0, right_spread=0.0, seed=seed,
                               time_budget=wall, polish=True))
        row["protocol_a"]["summary"] = summarise(row["protocol_a"]["runs"], exact)

        marker = (kind, tag[0], tag[1]) if kind == "synthetic" else None
        if marker in PROTOCOL_B:
            row["protocol_b"] = {"evaluation_budget": EVALUATION_BUDGET, "runs": []}
            for seed in SEEDS:
                row["protocol_b"]["runs"].append(
                    run_ga_matched(crisp(benchmark, threshold), 0.0, 0.0,
                                   left_spread=0.0, right_spread=0.0, seed=seed,
                                   evaluation_budget=EVALUATION_BUDGET, polish=True))
            row["protocol_b"]["summary"] = summarise(row["protocol_b"]["runs"], exact)

        rows.append(row)
        DESTINATION.write_text(json.dumps({
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "protocol_a": "GA wall-clock matched to the exact solve (floor 5 s), "
                          "restart-on-stagnation, equal-cost polish",
            "protocol_b": f"GA given {EVALUATION_BUDGET:,} evaluations, "
                          "restart-on-stagnation, equal-cost polish",
            "implementation_note": "The GA is pure Python and the exact solver is "
                                   "compiled C++, so wall-clock matching favours the "
                                   "exact solver; protocol B reports the "
                                   "evaluation-budget view alongside it.",
            "seeds": list(SEEDS),
            "rows": rows,
        }, indent=2))

        summary = row["protocol_a"]["summary"]
        print(f"{label:<20} exact={row['exact_cost']:>4} "
              f"A: cost {summary['cost_hits']}/5, "
              f"secondary {summary['secondary_hits_before_polish']}->"
              f"{summary['secondary_hits_after_polish']}/5 "
              f"(budget {wall:.1f}s, {summary['mean_evaluations']:.0f} evals)", flush=True)
        if "protocol_b" in row:
            s = row["protocol_b"]["summary"]
            print(f"{'':<20} B: cost {s['cost_hits']}/5, secondary "
                  f"{s['secondary_hits_before_polish']}->{s['secondary_hits_after_polish']}/5 "
                  f"({s['mean_evaluations']:.0f} evals, {s['mean_wall_seconds']:.0f}s)", flush=True)

    print(f"\nwritten: {DESTINATION}")


if __name__ == "__main__":
    main()
