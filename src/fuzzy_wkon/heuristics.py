"""Metaheuristic baseline for the finite-catalog design problem.

The fuzzy and interval allocation literature this study positions against
(Type-2 fuzzy, interval-dominance and fuzzy-clustering formulations) solves the
design problem with evolutionary or swarm search.  A referee is entitled to ask
what the exact certificate buys over that practice, so the comparison has to be
run rather than asserted.

The genetic algorithm below is deliberately conventional -- tournament
selection, uniform crossover, per-gene mutation, elitism, feasibility-first
comparison.  It is not tuned to lose: the population, generation count and
restart budget are scaled with the instance, and the same lexicographic
objective (minimum cost, then maximum reliability) is used, so the only
difference from the exact method is the absence of a certificate.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Any, Sequence

from .reliability import capped_distribution
from .solver import prepare_options


@dataclass
class GAResult:
    feasible: bool
    best_cost: int | None
    best_reliability: float | None
    best_chromosome: list[int] | None
    evaluations: int
    runtime_seconds: float
    generations_run: int


def _evaluate(chromosome: Sequence[int], rows, threshold: int) -> tuple[int, float]:
    selected: list[tuple[float, int]] = []
    cost = 0
    for row, choice_index in zip(rows, chromosome):
        option = row[choice_index]
        cost += option.cost
        if option.weight > 0 or option.reliability > 0.0:
            selected.append((option.reliability, option.weight))
    reliability = capped_distribution(selected, threshold)[threshold]
    return cost, reliability


@dataclass
class _Row:
    cost: int
    reliability: float
    weight: int


def run_ga(benchmark: dict[str, Any], alpha: float, gamma: float,
           left_spread: float = 0.1, right_spread: float = 0.1,
           population: int | None = None, generations: int | None = None,
           restarts: int = 3, seed: int = 0,
           time_budget: float | None = None,
           polish_equal_cost: bool = False) -> GAResult:
    from time import perf_counter

    started = perf_counter()
    options, threshold, _scale = prepare_options(benchmark, alpha, gamma,
                                                 left_spread, right_spread)
    # Index options by their catalog position so a gene is a catalog index.
    rows: list[list[_Row]] = []
    for row in options:
        ordered = sorted(row, key=lambda option: option.choice_index)
        rows.append([_Row(option.cost, option.reliability, option.weight_scaled)
                     for option in ordered])

    n = len(rows)
    target = benchmark["R_target"]
    population = population or max(60, 12 * n)
    generations = generations or max(200, 40 * n)

    rng = random.Random(seed)
    evaluations = 0
    generations_run = 0

    best: tuple[int, float, list[int]] | None = None

    def fitness_key(entry: tuple[int, float, list[int]]):
        cost, reliability, _ = entry
        feasible = reliability + 1e-12 >= target
        # Feasibility first, then minimum cost, then maximum reliability.
        return (0 if feasible else 1,
                cost if feasible else 0,
                -reliability if feasible else (target - reliability))

    def evaluate(chromosome: list[int]) -> tuple[int, float, list[int]]:
        nonlocal evaluations
        evaluations += 1
        cost, reliability = _evaluate(chromosome, rows, threshold)
        return (cost, reliability, chromosome)

    for restart in range(restarts):
        pool = [evaluate([rng.randrange(len(rows[i])) for i in range(n)])
                for _ in range(population)]
        pool.sort(key=fitness_key)
        stagnant = 0
        for generation in range(generations):
            generations_run += 1
            if time_budget is not None and perf_counter() - started > time_budget:
                break
            children: list[tuple[int, float, list[int]]] = [pool[0], pool[1]]  # elitism
            while len(children) < population:
                def tournament() -> list[int]:
                    contenders = [pool[rng.randrange(population)] for _ in range(3)]
                    contenders.sort(key=fitness_key)
                    return contenders[0][2]

                mother, father = tournament(), tournament()
                child = [mother[i] if rng.random() < 0.5 else father[i] for i in range(n)]
                for i in range(n):
                    if rng.random() < 1.5 / n:
                        child[i] = rng.randrange(len(rows[i]))
                children.append(evaluate(child))
            children.sort(key=fitness_key)
            if best is None or fitness_key(children[0]) < fitness_key(best):
                best = children[0]
                stagnant = 0
            else:
                stagnant += 1
            pool = children
            if stagnant > max(30, 5 * n):
                break
        if time_budget is not None and perf_counter() - started > time_budget:
            break

    if polish_equal_cost and best is not None:
        cost, reliability, chromosome = best
        if reliability + 1e-12 >= target:
            improved = True
            current_chrom = list(chromosome)
            current_rel = reliability
            target_cost = cost
            # 1-gene and 2-gene equal-cost polish
            while improved:
                improved = False
                for i in range(n):
                    orig_choice = current_chrom[i]
                    for choice in range(len(rows[i])):
                        if choice == orig_choice:
                            continue
                        candidate = list(current_chrom)
                        candidate[i] = choice
                        cand_cost, cand_rel, _ = evaluate(candidate)
                        if cand_cost == target_cost and cand_rel > current_rel + 1e-12:
                            current_chrom = candidate
                            current_rel = cand_rel
                            improved = True
                            break
                    if improved:
                        break
                if not improved:
                    for i in range(n):
                        for j in range(i + 1, n):
                            orig_i = current_chrom[i]
                            orig_j = current_chrom[j]
                            for ci in range(len(rows[i])):
                                for cj in range(len(rows[j])):
                                    if ci == orig_i and cj == orig_j:
                                        continue
                                    candidate = list(current_chrom)
                                    candidate[i] = ci
                                    candidate[j] = cj
                                    cand_cost, cand_rel, _ = evaluate(candidate)
                                    if cand_cost == target_cost and cand_rel > current_rel + 1e-12:
                                        current_chrom = candidate
                                        current_rel = cand_rel
                                        improved = True
                                        break
                                if improved:
                                    break
                            if improved:
                                break
            best = (target_cost, current_rel, current_chrom)

    runtime = perf_counter() - started
    if best is None:
        return GAResult(False, None, None, None, evaluations, runtime, generations_run)
    cost, reliability, chromosome = best
    if reliability + 1e-12 < target:
        return GAResult(False, None, None, None, evaluations, runtime, generations_run)
    return GAResult(True, cost, reliability, list(chromosome), evaluations, runtime,
                    generations_run)


# ---------------------------------------------------------------------------
# Budget-matched GA with equal-cost polish (round-2 reviewer response).
#
# Two reviewers objected that the round-1 comparison starved the heuristic:
# the exact solver ran in C++ for 371 s while the GA ran in Python for 15 s and
# exited on its own stagnation rule.  The functions below remove every part of
# that objection that can be removed without rewriting the GA in C++:
#
#   * the stagnation rule restarts the search instead of ending it, so the GA
#     uses its whole budget;
#   * the budget can be given as wall-clock (matched to the exact solve) or as
#     an evaluation count (which is the currency that does not penalise Python
#     for being slower per candidate than the compiled solver);
#   * an equal-cost local search polishes the incumbent, which is the fix a
#     practitioner would apply to a cost-driven search whose secondary
#     objective is reliability.
#
# The language asymmetry is stated, not hidden: wall-clock matching favours the
# exact solver, so the evaluation-budget run is reported alongside it.
# ---------------------------------------------------------------------------


def equal_cost_polish(chromosome: list[int], rows, threshold: int, target: float,
                      max_rounds: int = 200) -> tuple[list[int], int, float, int]:
    """Steepest-ascent local search over designs of identical total cost.

    Two neighbourhoods are explored: single-subsystem substitutions that leave
    the cost unchanged, and compensating pairwise swaps whose cost deltas
    cancel.  Only feasible, strictly more reliable designs are accepted, so the
    polish can never worsen either objective.
    """

    current = list(chromosome)
    cost, reliability = _evaluate(current, rows, threshold)
    evaluations = 0
    n = len(rows)

    for _ in range(max_rounds):
        best_gain = 0.0
        best_move = None

        # Neighbourhood 1: one subsystem, same cost.
        for i in range(n):
            for j, option in enumerate(rows[i]):
                if j == current[i] or option.cost != rows[i][current[i]].cost:
                    continue
                trial = list(current)
                trial[i] = j
                _, trial_r = _evaluate(trial, rows, threshold)
                evaluations += 1
                if trial_r >= target and trial_r - reliability > best_gain + 1e-15:
                    best_gain, best_move = trial_r - reliability, trial

        # Neighbourhood 2: two subsystems, cancelling cost deltas.
        for i in range(n):
            for j, oi in enumerate(rows[i]):
                if j == current[i]:
                    continue
                delta = oi.cost - rows[i][current[i]].cost
                for k in range(i + 1, n):
                    for l, ok in enumerate(rows[k]):
                        if l == current[k] or ok.cost - rows[k][current[k]].cost != -delta:
                            continue
                        trial = list(current)
                        trial[i], trial[k] = j, l
                        _, trial_r = _evaluate(trial, rows, threshold)
                        evaluations += 1
                        if trial_r >= target and trial_r - reliability > best_gain + 1e-15:
                            best_gain, best_move = trial_r - reliability, trial

        if best_move is None:
            break
        current = best_move
        cost, reliability = _evaluate(current, rows, threshold)

    return current, cost, reliability, evaluations


def run_ga_matched(benchmark: dict[str, Any], alpha: float, gamma: float,
                   left_spread: float = 0.1, right_spread: float = 0.1,
                   seed: int = 0, time_budget: float | None = None,
                   evaluation_budget: int | None = None,
                   polish: bool = True) -> dict[str, Any]:
    """GA that spends its whole declared budget, then polishes at equal cost."""

    from time import perf_counter

    started = perf_counter()
    options, threshold, _scale = prepare_options(benchmark, alpha, gamma,
                                                 left_spread, right_spread)
    rows: list[list[_Row]] = []
    for row in options:
        ordered = sorted(row, key=lambda option: option.choice_index)
        rows.append([_Row(option.cost, option.reliability, option.weight_scaled)
                     for option in ordered])

    n = len(rows)
    target = benchmark["R_target"]
    population = max(60, 12 * n)
    rng = random.Random(seed)
    evaluations = 0
    restarts = 0
    best: tuple[int, float, list[int]] | None = None

    def exhausted() -> bool:
        if time_budget is not None and perf_counter() - started >= time_budget:
            return True
        if evaluation_budget is not None and evaluations >= evaluation_budget:
            return True
        return False

    def key(entry):
        cost, reliability, _ = entry
        feasible = reliability + 1e-12 >= target
        return (0 if feasible else 1,
                cost if feasible else 0,
                -reliability if feasible else (target - reliability))

    def evaluate(chromosome: list[int]):
        nonlocal evaluations
        evaluations += 1
        cost, reliability = _evaluate(chromosome, rows, threshold)
        return (cost, reliability, chromosome)

    while not exhausted():
        restarts += 1
        pool = [evaluate([rng.randrange(len(rows[i])) for i in range(n)])
                for _ in range(population)]
        pool.sort(key=key)
        stagnant = 0
        while not exhausted() and stagnant <= max(30, 5 * n):
            children = [pool[0], pool[1]]
            while len(children) < population and not exhausted():
                def tournament() -> list[int]:
                    contenders = [pool[rng.randrange(len(pool))] for _ in range(3)]
                    contenders.sort(key=key)
                    return contenders[0][2]

                mother, father = tournament(), tournament()
                child = [mother[i] if rng.random() < 0.5 else father[i] for i in range(n)]
                for i in range(n):
                    if rng.random() < 1.5 / n:
                        child[i] = rng.randrange(len(rows[i]))
                children.append(evaluate(child))
            children.sort(key=key)
            if best is None or key(children[0]) < key(best):
                best, stagnant = children[0], 0
            else:
                stagnant += 1
            pool = children

    result: dict[str, Any] = {
        "seed": seed,
        "evaluations": evaluations,
        "restarts": restarts,
        "wall_seconds": perf_counter() - started,
        "time_budget": time_budget,
        "evaluation_budget": evaluation_budget,
    }
    if best is None or best[1] + 1e-12 < target:
        result.update({"feasible": False, "cost": None, "reliability": None,
                       "chromosome": None, "polished": False})
        return result

    cost, reliability, chromosome = best
    result.update({"feasible": True, "cost": cost,
                   "reliability_before_polish": reliability,
                   "chromosome_before_polish": list(chromosome)})

    if polish:
        chromosome, cost, reliability, polish_evals = equal_cost_polish(
            chromosome, rows, threshold, target)
        result["polish_evaluations"] = polish_evals
        result["polish_gain"] = reliability - result["reliability_before_polish"]
    result.update({"reliability": reliability, "chromosome": list(chromosome),
                   "polished": bool(polish)})
    return result
