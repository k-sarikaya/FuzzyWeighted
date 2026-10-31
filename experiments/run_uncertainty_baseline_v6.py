"""Oracle-substitution experiment: possibility theory vs Liu's uncertainty theory.

The reviewers of the submitted version asked what the possibilistic machinery buys
relative to a competing epistemic framework. The honest way to answer is to hold
the decision problem fixed -- same catalogs, same cost objective, same target --
and swap only the object that grades a design's reliability.

Arm P: the possibilistic evaluator of the main text, at the modal threshold.
Arm U: the uncertain-measure evaluator of Wei et al. (2025), Symmetry 17:908,
       with component reliabilities read as uncertain state measures and the
       capacity TFNs carried over as zigzag uncertain variables.

Both arms are solved to certified optimality by the same lexicographic branch and
bound (minimum cost, then maximum reliability), so any difference in the selected
design is attributable to the calculus and not to the search.

    PYTHONPATH=src python3 experiments/run_uncertainty_baseline_v6.py
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from math import prod
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS                       # noqa: E402
from fuzzy_wkon.uncertain_oracle import (                          # noqa: E402
    evaluate_design_uncertain, optimistic_bound_uncertain, verify_example)

DESTINATION = ROOT / "results" / "generated" / "uncertainty_baseline_v6.json"
SPREAD = 0.1
RM_TOLERANCE = 1e-6


def solve_uncertain(benchmark: dict, threshold: float, target: float) -> dict:
    """Lexicographic exact search under the uncertain-measure oracle."""
    rows = [sorted(range(len(r)), key=lambda j: r[j][1]) for r in benchmark["choices"]]
    n = len(rows)
    suffix_min_cost = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        suffix_min_cost[i] = suffix_min_cost[i + 1] + min(
            option[1] for option in benchmark["choices"][i])

    best = {"cost": None, "rm": None, "design": None}
    stats = {"nodes": 0, "leaves": 0, "pruned_cost": 0, "pruned_bound": 0}

    def recurse(i: int, prefix: list[int], cost: int) -> None:
        stats["nodes"] += 1
        if best["cost"] is not None and cost + suffix_min_cost[i] > best["cost"]:
            stats["pruned_cost"] += 1
            return
        if i == n:
            stats["leaves"] += 1
            _, rm = evaluate_design_uncertain(benchmark, prefix, threshold,
                                              SPREAD, SPREAD)
            if rm + RM_TOLERANCE < target:
                return
            if (best["cost"] is None or cost < best["cost"]
                    or (cost == best["cost"] and rm > best["rm"] + 1e-12)):
                best.update(cost=cost, rm=rm, design=list(prefix))
            return
        bound = optimistic_bound_uncertain(benchmark, prefix, threshold,
                                           SPREAD, SPREAD)
        if bound + RM_TOLERANCE < target:
            stats["pruned_bound"] += 1
            return
        for choice in rows[i]:
            recurse(i + 1, prefix + [choice], cost + benchmark["choices"][i][choice][1])

    started = perf_counter()
    sys.setrecursionlimit(10000)
    recurse(0, [], 0)
    return {**best, "stats": stats, "seconds": perf_counter() - started}


def solve_possibilistic(benchmark: dict, threshold: int) -> dict:
    """The main text's certified solver, run at the same crisp threshold."""
    from fuzzy_wkon.cpp_backend import solve_at_integer_threshold_cpp
    started = perf_counter()
    out = solve_at_integer_threshold_cpp(benchmark, threshold, time_limit=1200.0)
    return {
        "cost": out.get("best_cost"),
        "reliability": out.get("best_reliability"),
        "design": out.get("best_chromosome"),
        "status": out.get("status"),
        "seconds": perf_counter() - started,
    }


def main() -> None:
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "oracle substitution: possibility theory vs Liu's uncertainty theory",
        "calculus_note": ("The map p_i -> M{x_i = 1} and TFN -> zigzag uncertain "
                          "variable is a translation between formalisms, not an "
                          "identity; an uncertain measure and a possibility "
                          "distribution are different objects."),
        "spread": SPREAD,
        "oracle_verification": verify_example(),
        "rows": [],
    }

    for name in ("pharma_adapted", "fyffe_adapted"):
        benchmark = BENCHMARKS[name]
        threshold = float(benchmark["K_tfn"][1])          # modal threshold
        target = benchmark["R_target"]
        space = prod(len(r) for r in benchmark["choices"])

        arm_u = solve_uncertain(benchmark, threshold, target)
        arm_p = solve_possibilistic(benchmark, int(threshold))

        same = (arm_u["design"] is not None and arm_p["design"] is not None
                and list(arm_u["design"]) == list(arm_p["design"]))
        row = {
            "benchmark": name,
            "n": len(benchmark["choices"]),
            "search_space": space,
            "threshold": threshold,
            "target": target,
            "possibilistic": arm_p,
            "uncertain": arm_u,
            "same_design": same,
            "cost_gap": (None if arm_u["cost"] is None or arm_p["cost"] is None
                         else arm_u["cost"] - arm_p["cost"]),
        }
        report["rows"].append(row)
        print(f"{name:<16} space={space:.3e}  "
              f"P: cost={arm_p['cost']} R={arm_p['reliability']} ({arm_p['seconds']:.2f}s)  "
              f"U: cost={arm_u['cost']} RM={arm_u['rm']} ({arm_u['seconds']:.1f}s, "
              f"{arm_u['stats']['leaves']} leaves)  same_design={same}", flush=True)
        DESTINATION.write_text(json.dumps(report, indent=2, default=str))

    print(f"\nwritten: {DESTINATION}")


if __name__ == "__main__":
    main()
