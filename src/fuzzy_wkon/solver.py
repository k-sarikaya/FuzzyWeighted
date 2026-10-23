"""Two-phase exact solver for the V2 objective.

Objective: minimum feasible cost, then maximum reliability among designs at
that minimum cost.  Zimmermann satisfaction is deliberately not used by the
solver objective.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import lcm
from time import perf_counter
from typing import Any

from .fuzzy import TFN
from .reliability import capped_distribution


@dataclass(frozen=True)
class Option:
    choice_index: int
    reliability: float
    cost: int
    weight_scaled: int


def _fuzzy_weight(weight: float, alpha: float, gamma: float, left_spread: float, right_spread: float) -> Fraction:
    w = Fraction(str(weight))
    a = Fraction(str(alpha))
    g = Fraction(str(gamma))
    left = Fraction(str(left_spread))
    right = Fraction(str(right_spread))
    lower = w * (1 - left) + a * left * w
    upper = w * (1 + right) - a * right * w
    return lower + g * (upper - lower)


def _fuzzy_threshold(tfn_values: tuple[float, float, float], alpha: float, gamma: float) -> Fraction:
    a = Fraction(str(alpha))
    g = Fraction(str(gamma))
    k1, k2, k3 = (Fraction(str(value)) for value in tfn_values)
    lower = k1 + a * (k2 - k1)
    upper = k3 - a * (k3 - k2)
    return upper - g * (upper - lower)


def prepare_options(benchmark: dict[str, Any], alpha: float, gamma: float,
                    left_spread: float = 0.1, right_spread: float = 0.1,
                    unit_spreads: list[list[tuple[float, float]]] | None = None) -> tuple[list[list[Option]], int, int]:
    threshold = _fuzzy_threshold(benchmark["K_tfn"], alpha, gamma)
    fractions: list[tuple[list[tuple[int, float, int, Fraction]], ...]] = []
    denominators = [threshold.denominator]
    for i, choices in enumerate(benchmark["choices"]):
        row: list[tuple[int, float, int, Fraction]] = []
        for index, (reliability, cost, weight) in enumerate(choices):
            if index == 0:
                weight_fraction = Fraction(0, 1)
            else:
                if unit_spreads is not None and i < len(unit_spreads) and index < len(unit_spreads[i]):
                    opt_l, opt_r = unit_spreads[i][index]
                else:
                    opt_l, opt_r = left_spread, right_spread
                weight_fraction = _fuzzy_weight(weight, alpha, gamma, opt_l, opt_r)
            denominators.append(weight_fraction.denominator)
            row.append((index, float(reliability), int(cost), weight_fraction))
        fractions.append(tuple(row))
    scale = 1
    for denominator in denominators:
        scale = lcm(scale, denominator)
    processed: list[list[Option]] = []
    for row in fractions:
        options = [Option(i, p, c, int(weight * scale)) for i, p, c, weight in row]
        options.sort(key=lambda option: (option.cost, -option.reliability, -option.weight_scaled, option.choice_index))
        processed.append(options)
    return processed, int(threshold * scale), scale


def evaluate_design(chromosome: list[int], benchmark: dict[str, Any], alpha: float, gamma: float,
                    left_spread: float = 0.1, right_spread: float = 0.1) -> dict[str, Any]:
    options, threshold, scale = prepare_options(benchmark, alpha, gamma, left_spread, right_spread)
    selected: list[tuple[float, int]] = []
    cost = 0
    for row, choice_index in zip(options, chromosome):
        option = next(option for option in row if option.choice_index == choice_index)
        cost += option.cost
        if option.choice_index != 0:
            selected.append((option.reliability, option.weight_scaled))
    pmf = capped_distribution(selected, threshold)
    return {
        "chromosome": list(chromosome),
        "cost": cost,
        "reliability": pmf[-1],
        "feasible": pmf[-1] + 1e-12 >= benchmark["R_target"],
        "threshold_scaled": threshold,
        "scale": scale,
    }


def _suffix_min_cost(options: list[list[Option]], threshold: int) -> list[list[int]]:
    inf = 10**12
    suffix = [[inf] * (threshold + 1) for _ in range(len(options) + 1)]
    suffix[-1][0] = 0
    for index in range(len(options) - 1, -1, -1):
        for need in range(threshold + 1):
            suffix[index][need] = min(
                option.cost + suffix[index + 1][max(0, need - option.weight_scaled)]
                for option in options[index]
            )
    return suffix


def _suffix_max_weight(options: list[list[Option]]) -> list[int]:
    suffix = [0] * (len(options) + 1)
    for index in range(len(options) - 1, -1, -1):
        suffix[index] = suffix[index + 1] + max(option.weight_scaled for option in options[index])
    return suffix


def _optimistic_reliability(pmf: dict[int, float], options: list[list[Option]], start: int, threshold: int) -> float:
    bound = pmf.copy()
    for row in options[start:]:
        max_p = max(option.reliability for option in row)
        max_w = max(option.weight_scaled for option in row)
        bound = capped_distribution_from_pmf(bound, max_p, max_w, threshold)
    return bound.get(threshold, 0.0)


def capped_distribution_from_pmf(pmf: dict[int, float], probability: float, weight: int, threshold: int) -> dict[int, float]:
    nxt: dict[int, float] = {}
    for current, mass in pmf.items():
        nxt[current] = nxt.get(current, 0.0) + mass * (1.0 - probability)
        target = min(threshold, current + weight)
        nxt[target] = nxt.get(target, 0.0) + mass * probability
    return nxt


def _heuristic_incumbent(options: list[list[Option]], threshold: int, target: float,
                         beam_width: int = 200) -> tuple[int, list[int]] | None:
    """Find a feasible cost upper bound without affecting exactness."""

    states: list[tuple[int, dict[int, float], list[int], float]] = [(0, {0: 1.0}, [], 1.0)]
    for index, row in enumerate(options):
        candidates: list[tuple[int, dict[int, float], list[int], float]] = []
        for cost, pmf, chromosome, _ in states:
            for option in row:
                nxt = capped_distribution_from_pmf(pmf, option.reliability, option.weight_scaled, threshold)
                upper = _optimistic_reliability(nxt, options, index + 1, threshold)
                if upper + 1e-12 < target:
                    continue
                candidates.append((cost + option.cost, nxt, chromosome + [option.choice_index], upper))
        candidates.sort(key=lambda state: (state[0], -state[3]))
        states = candidates[:beam_width]
        if not states:
            return None
    feasible = [state for state in states if state[1].get(threshold, 0.0) + 1e-12 >= target]
    if not feasible:
        return None
    best = min(feasible, key=lambda state: (state[0], -state[1].get(threshold, 0.0)))
    return best[0], best[2]


def solve_exact(benchmark: dict[str, Any], alpha: float, gamma: float,
                left_spread: float = 0.1, right_spread: float = 0.1) -> dict[str, Any]:
    options, threshold, scale = prepare_options(benchmark, alpha, gamma, left_spread, right_spread)
    suffix = _suffix_min_cost(options, threshold)
    max_weight_suffix = _suffix_max_weight(options)
    initial = {0: 1.0}
    stats = {"nodes": 0, "leaves": 0, "pruned_cost": 0, "pruned_capacity": 0, "pruned_reliability": 0}
    start_clock = perf_counter()

    def selection_state(chromosome: list[int]) -> tuple[int, dict[int, float]]:
        cost = 0
        pmf = {0: 1.0}
        for row, choice_index in zip(options, chromosome):
            option = next(item for item in row if item.choice_index == choice_index)
            cost += option.cost
            pmf = capped_distribution_from_pmf(pmf, option.reliability, option.weight_scaled, threshold)
        return cost, pmf

    # A feasible incumbent is only a pruning aid; optimality still comes from
    # the complete admissible search below. The deterministic beam heuristic
    # supplies a useful upper bound for larger catalogs such as Fyffe.
    beam_width = 5000 if len(options) > 10 else 200
    heuristic = _heuristic_incumbent(options, threshold, benchmark["R_target"], beam_width)
    best_cost = heuristic[0] if heuristic is not None else float("inf")

    def phase_one(index: int, cost: int, capacity: int, pmf: dict[int, float], chromosome: list[int]) -> None:
        nonlocal best_cost
        stats["nodes"] += 1
        if cost > best_cost:
            stats["pruned_cost"] += 1
            return
        if capacity + max_weight_suffix[index] < threshold:
            stats["pruned_capacity"] += 1
            return
        if _optimistic_reliability(pmf, options, index, threshold) + 1e-12 < benchmark["R_target"]:
            stats["pruned_reliability"] += 1
            return
        # suffix cost is a valid lower bound even though this first pass
        # deliberately keeps the capacity component conservative.
        lb = suffix[index][max(0, threshold - capacity)]
        if cost + lb > best_cost:
            stats["pruned_capacity"] += 1
            return
        if index == len(options):
            stats["leaves"] += 1
            if pmf.get(threshold, 0.0) + 1e-12 >= benchmark["R_target"]:
                best_cost = min(best_cost, cost)
            return
        for option in options[index]:
            next_cost = cost + option.cost
            if next_cost > best_cost:
                stats["pruned_cost"] += 1
                continue
            nxt = capped_distribution_from_pmf(pmf, option.reliability, option.weight_scaled, threshold)
            chromosome.append(option.choice_index)
            phase_one(index + 1, next_cost, capacity + option.weight_scaled, nxt, chromosome)
            chromosome.pop()

    # The phase-one capacity LB is intentionally supplemented by exhaustive
    # feasibility traversal; reliability cannot be safely inferred from weight
    # alone.  This remains exact and keeps the proof simple.
    phase_one(0, 0, 0, initial, [])
    if best_cost == float("inf"):
        return {
            "feasible": False,
            "best_cost": None,
            "best_reliability": None,
            "best_chromosome": None,
            "scale": scale,
            "threshold_scaled": threshold,
            "runtime_seconds": perf_counter() - start_clock,
            "stats": stats,
        }

    best_reliability = -1.0
    best_chromosome: list[int] | None = None

    def phase_two(index: int, cost: int, capacity: int, pmf: dict[int, float], chromosome: list[int]) -> None:
        nonlocal best_reliability, best_chromosome
        stats["nodes"] += 1
        if cost > best_cost:
            stats["pruned_cost"] += 1
            return
        if capacity + max_weight_suffix[index] < threshold:
            stats["pruned_capacity"] += 1
            return
        if cost + suffix[index][max(0, threshold - capacity)] > best_cost:
            stats["pruned_capacity"] += 1
            return
        upper = _optimistic_reliability(pmf, options, index, threshold)
        if upper + 1e-12 < benchmark["R_target"]:
            stats["pruned_reliability"] += 1
            return
        if index == len(options):
            stats["leaves"] += 1
            reliability = pmf.get(threshold, 0.0)
            if cost == best_cost and reliability + 1e-12 >= benchmark["R_target"]:
                if reliability > best_reliability + 1e-12:
                    best_reliability = reliability
                    best_chromosome = chromosome.copy()
            return
        for option in options[index]:
            next_cost = cost + option.cost
            if next_cost > best_cost:
                stats["pruned_cost"] += 1
                continue
            nxt = capped_distribution_from_pmf(pmf, option.reliability, option.weight_scaled, threshold)
            chromosome.append(option.choice_index)
            phase_two(index + 1, next_cost, capacity + option.weight_scaled, nxt, chromosome)
            chromosome.pop()

    phase_two(0, 0, 0, initial, [])
    elapsed = perf_counter() - start_clock
    return {
        "feasible": True,
        "best_cost": int(best_cost),
        "best_reliability": float(best_reliability),
        "best_chromosome": best_chromosome,
        "alpha": alpha,
        "gamma": gamma,
        "left_spread": left_spread,
        "right_spread": right_spread,
        "scale": scale,
        "threshold_scaled": threshold,
        "runtime_seconds": elapsed,
        "objective": ["min_cost", "max_reliability"],
        "stats": stats,
    }
