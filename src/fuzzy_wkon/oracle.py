"""Independent exhaustive oracle for small adapted instances."""

from __future__ import annotations

from itertools import product

import numpy as np

from .reliability import capped_distribution
from .solver import prepare_options
from .solver import capped_distribution_from_pmf


def exhaustive_optimum(benchmark: dict, alpha: float, gamma: float,
                       left_spread: float = 0.1, right_spread: float = 0.1) -> dict:
    best: dict | None = None
    options, threshold, scale = prepare_options(benchmark, alpha, gamma, left_spread, right_spread)
    for choices in product(*options):
        cost = sum(option.cost for option in choices)
        selected = [(option.reliability, option.weight_scaled) for option in choices if option.choice_index != 0]
        reliability = capped_distribution(selected, threshold)[-1]
        if reliability + 1e-12 < benchmark["R_target"]:
            continue
        result = {
            "feasible": True,
            "cost": cost,
            "reliability": reliability,
            "chromosome": [option.choice_index for option in choices],
            "threshold_scaled": threshold,
            "scale": scale,
        }
        key = (cost, -reliability)
        if best is None or key < (best["cost"], -best["reliability"]):
            best = result
    return best or {"feasible": False, "cost": None, "reliability": None, "chromosome": None}


def vectorized_exhaustive_optimum(benchmark: dict, alpha: float, gamma: float,
                                  left_spread: float = 0.1, right_spread: float = 0.1,
                                  state_chunk: int = 64) -> dict:
    """Independent exhaustive oracle using batched state enumeration.

    Every catalog chromosome is still evaluated. NumPy only batches the
    independent component states so this remains an oracle rather than a
    branch-and-bound implementation.
    """

    options, threshold, scale = prepare_options(benchmark, alpha, gamma, left_spread, right_spread)
    grids = np.meshgrid(*[np.arange(len(row), dtype=np.int16) for row in options], indexing="ij")
    choice_matrix = np.stack([grid.reshape(-1) for grid in grids], axis=1)
    count, n = choice_matrix.shape
    cost_matrix = np.zeros((count, n), dtype=np.float64)
    weight_matrix = np.zeros((count, n), dtype=np.int64)
    probability_matrix = np.zeros((count, n), dtype=np.float64)
    for component, row in enumerate(options):
        cost_values = np.asarray([item.cost for item in row], dtype=np.float64)
        weight_values = np.asarray([item.weight_scaled for item in row], dtype=np.int64)
        probability_values = np.asarray([item.reliability for item in row], dtype=np.float64)
        indexes = choice_matrix[:, component]
        cost_matrix[:, component] = cost_values[indexes]
        weight_matrix[:, component] = weight_values[indexes]
        probability_matrix[:, component] = probability_values[indexes]
    costs = cost_matrix.sum(axis=1)
    reliabilities = np.zeros(count, dtype=np.float64)
    for start in range(0, 1 << n, state_chunk):
        stop = min(1 << n, start + state_chunk)
        masks = np.arange(start, stop, dtype=np.int16)
        states = ((masks[:, None] >> np.arange(n, dtype=np.int16)) & 1).astype(bool)
        state_weights = weight_matrix @ states.T
        success = state_weights >= threshold
        p = probability_matrix[:, None, :]
        probabilities = np.where(states[None, :, :], p, 1.0 - p).prod(axis=2)
        reliabilities += (probabilities * success).sum(axis=1)
    feasible = reliabilities + 1e-12 >= benchmark["R_target"]
    if not np.any(feasible):
        return {"feasible": False, "cost": None, "reliability": None, "chromosome": None}
    feasible_indexes = np.flatnonzero(feasible)
    best_cost = costs[feasible_indexes].min()
    cost_indexes = feasible_indexes[costs[feasible_indexes] == best_cost]
    best_index = cost_indexes[np.argmax(reliabilities[cost_indexes])]
    return {
        "feasible": True,
        "cost": int(best_cost),
        "reliability": float(reliabilities[best_index]),
        "chromosome": choice_matrix[best_index].astype(int).tolist(),
        "threshold_scaled": threshold,
        "scale": scale,
    }


def pareto_dp_optimum(benchmark: dict, alpha: float, gamma: float,
                      left_spread: float = 0.1, right_spread: float = 0.1,
                      max_cost: int = 45) -> dict:
    """Exact cost-bounded DP with first-order stochastic dominance pruning.

    At a fixed subsystem index and accumulated cost, a PMF whose tail
    probabilities dominate another PMF remains dominant after every common
    non-negative future convolution. Therefore dominated PMFs can be removed
    without changing the best final reliability.
    """

    options, threshold, scale = prepare_options(benchmark, alpha, gamma, left_spread, right_spread)

    def dominates(first: dict[int, float], second: dict[int, float]) -> bool:
        tail_first = 0.0
        tail_second = 0.0
        for weight in range(threshold, -1, -1):
            tail_first += first.get(weight, 0.0)
            tail_second += second.get(weight, 0.0)
            if tail_first + 1e-12 < tail_second:
                return False
        return True

    states: dict[int, list[tuple[dict[int, float], list[int]]]] = {0: [({0: 1.0}, [])]}
    for index, row in enumerate(options):
        expanded: dict[int, list[tuple[dict[int, float], list[int]]]] = {}
        for cost, items in states.items():
            for pmf, chromosome in items:
                for option in row:
                    next_cost = cost + option.cost
                    if next_cost > max_cost:
                        continue
                    next_pmf = capped_distribution_from_pmf(pmf, option.reliability, option.weight_scaled, threshold)
                    expanded.setdefault(next_cost, []).append((next_pmf, chromosome + [option.choice_index]))
        states = {}
        for cost, items in expanded.items():
            kept: list[tuple[dict[int, float], list[int]]] = []
            for pmf, chromosome in items:
                if any(dominates(existing, pmf) for existing, _ in kept):
                    continue
                kept = [(existing, old_chromosome) for existing, old_chromosome in kept if not dominates(pmf, existing)]
                kept.append((pmf, chromosome))
            states[cost] = kept

    best: dict | None = None
    for cost, items in states.items():
        for pmf, chromosome in items:
            reliability = pmf.get(threshold, 0.0)
            if reliability + 1e-12 < benchmark["R_target"]:
                continue
            candidate = {"feasible": True, "cost": cost, "reliability": reliability,
                         "chromosome": chromosome, "threshold_scaled": threshold, "scale": scale,
                         "engine": "pareto_cost_pmf_dp"}
            if best is None or (candidate["cost"], -candidate["reliability"]) < (best["cost"], -best["reliability"]):
                best = candidate
    return best or {"feasible": False, "cost": None, "reliability": None, "chromosome": None,
                    "threshold_scaled": threshold, "scale": scale, "engine": "pareto_cost_pmf_dp"}
