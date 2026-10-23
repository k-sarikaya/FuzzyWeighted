"""Reliability oracle of Wei et al. (2025), Symmetry 17:908, under Liu's uncertainty theory.

This module exists so that the possibilistic framework of the main text can be
compared against a *different epistemic calculus* on identical instances, rather
than against a rhetorical description of one. Wei et al. evaluate a weighted
k-out-of-n system in which the component weights and states are independent
uncertain variables and the threshold is a crisp constant. Their Theorem 4 turns
the reliability measure into the root of a one-dimensional monotone condition,

    RM = sup { nu in [0,1] : sum_i Phi_i^{-1}(1 - nu) * 1_i(nu) >= k },

where ``Phi_i`` is the uncertainty distribution of the weight of component ``i``
and ``1_i(nu)`` indicates whether that component is counted at confidence ``nu``.
For the binary case ``1_i(nu) = 1`` exactly when ``nu <= v_i``, with
``v_i = M{x_i = 1}`` the uncertain reliability measure of the component state.
(The published statement of the step function is written the other way round;
the multi-state companion paper, Symmetry 17:912, states the same object with the
orientation used here, and only this orientation reproduces the worked example.)

Two uncertainty distributions are supported. ``linear`` is Liu's L(a, b), which
Wei et al. use in their example. ``zigzag`` is Liu's Z(a, b, c), which carries the
same three parameters as a triangular fuzzy number and is therefore the natural
image of the capacity TFNs of the main text under a change of calculus. The
mapping is a translation between formalisms, not an identity: an uncertain
measure and a possibility distribution are different objects, and nothing here
asserts otherwise.

Verification against the published example lives in ``verify_example()``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

__all__ = ["UncertainWeight", "reliability_measure", "verify_example"]


@dataclass(frozen=True)
class UncertainWeight:
    """An uncertain weight with an invertible distribution.

    ``linear``: L(a, b)      -> Phi^{-1}(p) = a + p (b - a)
    ``zigzag``: Z(a, b, c)   -> Phi^{-1}(p) = a + 2p(b - a)        for p < 1/2
                                            = 2b - c + 2p(c - b)   for p >= 1/2
    """

    kind: Literal["linear", "zigzag"]
    params: tuple[float, ...]

    def __post_init__(self) -> None:
        need = 2 if self.kind == "linear" else 3
        if len(self.params) != need:
            raise ValueError(f"{self.kind} weight needs {need} parameters")
        if list(self.params) != sorted(self.params):
            raise ValueError("distribution parameters must be non-decreasing")

    def inverse(self, p: float) -> float:
        if not 0.0 <= p <= 1.0:
            raise ValueError("inverse distribution takes a probability level")
        if self.kind == "linear":
            a, b = self.params
            return a + p * (b - a)
        a, b, c = self.params
        return a + 2.0 * p * (b - a) if p < 0.5 else 2.0 * b - c + 2.0 * p * (c - b)


def _attained(nu: float, weights: Sequence[UncertainWeight],
              states: Sequence[float]) -> float:
    """Total weight credited at confidence ``nu``: monotone non-increasing in ``nu``."""
    level = 1.0 - nu
    return sum(w.inverse(level) for w, v in zip(weights, states) if nu <= v)


def reliability_measure(weights: Sequence[UncertainWeight], states: Sequence[float],
                        threshold: float, tolerance: float = 1e-9) -> float:
    """Return ``RM`` by bisection on the condition of Theorem 4.

    ``_attained`` is non-increasing in ``nu`` (the inverse distribution falls and
    components drop out as ``nu`` rises), so the condition is a downward-closed
    set and bisection is exact up to ``tolerance``.
    """
    if len(weights) != len(states):
        raise ValueError("weights and states must have equal length")
    if _attained(1.0, weights, states) >= threshold:
        return 1.0
    if _attained(0.0, weights, states) < threshold:
        return 0.0
    low, high = 0.0, 1.0            # attained(low) >= k > attained(high)
    while high - low > tolerance:
        mid = 0.5 * (low + high)
        if _attained(mid, weights, states) >= threshold:
            low = mid
        else:
            high = mid
    return low


def verify_example() -> dict[str, object]:
    """Reproduce Example 2 of Wei et al. (2025), Symmetry 17:908.

    Returns the computed partial sums alongside the published ones. The published
    value at nu = 0.9 is reproduced exactly; the remaining three published sums
    are not reproducible from the paper's own Tables 3 and 4 by ANY subset of
    components (verified by exhaustive enumeration over all 2^10 subsets), so we
    treat them as arithmetic slips in the source. The reported reliability
    measure RM = 0.6 is reproduced.
    """
    states = [0.5, 0.6, 0.7, 0.5, 0.9, 0.9, 0.8, 0.9, 0.7, 0.8]
    weights = [UncertainWeight("linear", p) for p in
               [(1, 2.2), (2, 2.5), (2, 2.8), (2, 3), (2, 3.4),
                (3, 3.6), (3, 4), (3, 4), (3, 4), (3, 4.6)]]
    published_sums = {0.9: 8.3, 0.8: 15.4, 0.7: 21.64, 0.6: 24.72}
    computed = {nu: round(_attained(nu, weights, states), 4) for nu in published_sums}
    return {
        "published_sums": published_sums,
        "computed_sums": computed,
        "published_RM_traverse": 0.6,
        "published_RM_bisection": 0.5977,
        "computed_RM": reliability_measure(weights, states, 24.0),
        "traverse_RM_at_step_0.1": max(
            (nu for nu in [round(1 - 0.1 * i, 1) for i in range(11)]
             if _attained(nu, weights, states) >= 24.0), default=0.0),
    }


if __name__ == "__main__":
    import json
    print(json.dumps(verify_example(), indent=2, default=str))


# ---------------------------------------------------------------------------
# Design-level interface: the same catalog problem, evaluated under the
# uncertainty calculus instead of the possibilistic one.
# ---------------------------------------------------------------------------

def design_weights_and_states(benchmark: dict, chromosome: Sequence[int],
                              left_spread: float, right_spread: float
                              ) -> tuple[list[UncertainWeight], list[float]]:
    """Translate one design of the catalog problem into the uncertainty calculus.

    Each subsystem contributes the option it selected. The option's reliability
    ``p`` becomes the uncertain reliability measure of the component state, and
    its capacity TFN ``(w(1-dl), w, w(1+dr))`` becomes the zigzag uncertain
    variable with the same three parameters. Options carrying no capacity and no
    reliability are dropped, exactly as in the possibilistic evaluator.
    """
    weights: list[UncertainWeight] = []
    states: list[float] = []
    for row, choice in zip(benchmark["choices"], chromosome):
        p, _cost, w = row[choice]
        if w <= 0 and p <= 0.0:
            continue
        weights.append(UncertainWeight(
            "zigzag", (w * (1.0 - left_spread), float(w), w * (1.0 + right_spread))))
        states.append(float(p))
    return weights, states


def evaluate_design_uncertain(benchmark: dict, chromosome: Sequence[int],
                              threshold: float, left_spread: float = 0.1,
                              right_spread: float = 0.1) -> tuple[int, float]:
    """Return ``(cost, RM)`` for a complete design under the uncertainty calculus."""
    cost = sum(row[c][1] for row, c in zip(benchmark["choices"], chromosome))
    weights, states = design_weights_and_states(
        benchmark, chromosome, left_spread, right_spread)
    return cost, reliability_measure(weights, states, threshold)


def optimistic_bound_uncertain(benchmark: dict, prefix: Sequence[int],
                               threshold: float, left_spread: float = 0.1,
                               right_spread: float = 0.1) -> float:
    """Admissible upper bound on ``RM`` over every completion of ``prefix``.

    Each unassigned subsystem is replaced by a virtual option carrying the largest
    state measure and the largest capacity available in its catalog. Because the
    capacities within a subsystem are proportional TFNs sharing the spreads, the
    option of largest modal capacity dominates every other option pointwise in the
    inverse distribution; together with the monotonicity of the uncertain measure
    in each argument (Liu's Theorem 1), the substituted system attains at least as
    much weight at every confidence level, so its ``RM`` bounds every completion.
    """
    weights, states = design_weights_and_states(
        benchmark, list(prefix), left_spread, right_spread)
    for row in benchmark["choices"][len(prefix):]:
        best_p = max(option[0] for option in row)
        best_w = max(option[2] for option in row)
        if best_w <= 0 and best_p <= 0.0:
            continue
        weights.append(UncertainWeight(
            "zigzag", (best_w * (1.0 - left_spread), float(best_w),
                       best_w * (1.0 + right_spread))))
        states.append(float(best_p))
    return reliability_measure(weights, states, threshold)
