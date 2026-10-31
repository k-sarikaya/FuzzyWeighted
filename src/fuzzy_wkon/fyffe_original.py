"""Exact reproduction of the original Fyffe--Hines--Lee numerical problem.

This module deliberately does *not* use the V2 weighted k-out-of-n:G model.
It preserves the source paper's series architecture, per-unit reliabilities,
costs and weights, and its integer parallel-redundancy decision.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path


SOURCE_TABLE = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "original"
    / "fyffe_1968_table_i.csv"
)


@dataclass(frozen=True)
class FyffeAlternative:
    unit: int
    alternative: int
    reliability: float
    cost: int
    weight: int


def load_table_i(path: Path = SOURCE_TABLE) -> list[list[FyffeAlternative]]:
    """Load the faithful transcription of Fyffe et al.'s Table I."""
    rows: dict[int, list[FyffeAlternative]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            item = FyffeAlternative(
                unit=int(row["functional_unit"]),
                alternative=int(row["design_alternative"]),
                reliability=float(row["reliability_R"]),
                cost=int(row["unit_cost_C"]),
                weight=int(row["unit_weight_W"]),
            )
            rows.setdefault(item.unit, []).append(item)
    if sorted(rows) != list(range(1, 15)):
        raise ValueError("Fyffe Table I must contain exactly functional units 1--14")
    return [rows[unit] for unit in range(1, 15)]


def parallel_reliability(unit_reliability: float, count: int) -> float:
    """Reliability of ``count`` identical independent parallel units."""
    if count < 1:
        raise ValueError("A Fyffe functional unit requires at least one component")
    return 1.0 - (1.0 - unit_reliability) ** count


def solve_original_fyffe(cost_limit: int = 130, weight_limit: int = 170) -> dict:
    """Exactly solve the source paper's constrained series allocation problem.

    Dynamic-programming states are exact cost--weight pairs. At a given pair,
    retaining only the largest log reliability is lossless: every remaining
    unit has the same feasible future choices, and the final series
    reliability is the product of per-unit reliabilities.
    """
    if cost_limit < 1 or weight_limit < 1:
        raise ValueError("Cost and weight limits must be positive")

    states: dict[tuple[int, int], tuple[float, tuple[tuple[int, int], ...]]] = {
        (0, 0): (0.0, ())
    }
    state_counts: list[int] = []
    for alternatives in load_table_i():
        next_states: dict[
            tuple[int, int], tuple[float, tuple[tuple[int, int], ...]]
        ] = {}
        for alternative in alternatives:
            maximum_count = min(
                cost_limit // alternative.cost, weight_limit // alternative.weight
            )
            for count in range(1, maximum_count + 1):
                added_cost = count * alternative.cost
                added_weight = count * alternative.weight
                log_reliability = math.log(
                    parallel_reliability(alternative.reliability, count)
                )
                for (cost, weight), (score, design) in states.items():
                    key = (cost + added_cost, weight + added_weight)
                    if key[0] > cost_limit or key[1] > weight_limit:
                        continue
                    candidate = (score + log_reliability, design + ((alternative.alternative, count),))
                    incumbent = next_states.get(key)
                    if incumbent is None or candidate[0] > incumbent[0] + 1e-15:
                        next_states[key] = candidate
        states = next_states
        state_counts.append(len(states))

    if not states:
        raise RuntimeError("No feasible Fyffe allocation under the supplied limits")
    (best_cost, best_weight), (best_log_reliability, best_design) = max(
        states.items(), key=lambda item: item[1][0]
    )
    return {
        "cost_limit": cost_limit,
        "weight_limit": weight_limit,
        "best_cost": best_cost,
        "best_weight": best_weight,
        "best_reliability": math.exp(best_log_reliability),
        "design": [
            {"functional_unit": index + 1, "alternative": alternative, "parallel_units": count}
            for index, (alternative, count) in enumerate(best_design)
        ],
        "final_state_count": len(states),
        "state_counts_by_layer": state_counts,
        "method": "exact cost-weight dynamic program maximizing log series reliability",
    }
