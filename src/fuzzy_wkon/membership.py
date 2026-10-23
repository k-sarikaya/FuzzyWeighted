"""Possibility, necessity and Hurwicz membership degrees for the success event.

V2 declared a single Hurwicz scalarisation and deliberately removed every
possibility/necessity statement because the underlying ordering had never been
constructed.  This module supplies that missing construction, so the Hurwicz
rule becomes one *declared member* of a decision-rule family rather than the
only available semantics.

Family
------
For a subset of operating components with total fuzzy capacity ``W`` and a
fuzzy threshold ``K`` the success event receives a degree in [0, 1]:

``pos``   optimistic bound, ``Pos(W >= K)``
``nec``   pessimistic bound, ``Nec(W >= K)``
``hurwicz``  the V2 indicator ``1{W_H >= K_H}`` at a declared ``(alpha, gamma)``

All three are non-decreasing in the modal capacity, which is exactly the
property used by the admissible-bound and stochastic-dominance proofs, so the
existing exact search machinery carries over unchanged.

Closed forms
------------
For triangular ``A=(a1,a2,a3)`` and ``B=(b1,b2,b3)`` the standard sup-min
construction gives

``Pos(A >= B) = 1``                                       if ``a2 >= b2``
``Pos(A >= B) = (a3-b1) / ((a3-a2) + (b2-b1))``           otherwise, clipped

and, through ``Nec(A >= B) = 1 - Pos(A <= B)``,

``Nec(A >= B) = 0``                                       if ``a2 <= b2``
``Nec(A >= B) = 1 - (b3-a1) / ((b3-b2) + (a2-a1))``       otherwise, clipped.

These expressions are verified against a brute-force sup-min grid in
``tests/test_membership.py``; they are not asserted from the literature alone.
"""

from __future__ import annotations

from typing import Callable

from .fuzzy import TFN, alpha_cut

__all__ = [
    "possibility_geq",
    "necessity_geq",
    "credibility_geq",
    "hurwicz_indicator",
    "proportional_capacity_tfn",
    "membership_profile",
    "MEASURES",
]


def _clip_unit(value: float) -> float:
    return 0.0 if value < 0.0 else (1.0 if value > 1.0 else value)


def possibility_geq(a: TFN, b: TFN) -> float:
    """Return ``Pos(a >= b)`` for two triangular fuzzy numbers."""

    if a.mode >= b.mode:
        return 1.0
    denominator = (a.right - a.mode) + (b.mode - b.left)
    if denominator <= 0.0:
        # Both branches degenerate: the comparison is crisp at the modes.
        return 1.0 if a.right >= b.left else 0.0
    return _clip_unit((a.right - b.left) / denominator)


def necessity_geq(a: TFN, b: TFN) -> float:
    """Return ``Nec(a >= b) = 1 - Pos(a <= b)`` for two triangular fuzzy numbers."""

    if a.mode <= b.mode:
        return 0.0
    denominator = (b.right - b.mode) + (a.mode - a.left)
    if denominator <= 0.0:
        return 1.0 if a.left >= b.right else 0.0
    return _clip_unit(1.0 - (b.right - a.left) / denominator)


def credibility_geq(a: TFN, b: TFN) -> float:
    """Liu's self-dual credibility measure, ``Cr = (Pos + Nec) / 2``.

    Cr is self-dual and, for this comparison, ``Cr >= 1/2`` reduces to the
    modal comparison ``a.mode >= b.mode``.  Note that the central Hurwicz
    column reproduces that test only when the capacity spreads are symmetric;
    see the corollary and its audit.
    """

    return 0.5 * (possibility_geq(a, b) + necessity_geq(a, b))


def hurwicz_indicator(a: TFN, b: TFN, alpha: float, gamma: float) -> float:
    """Return the V2 crisp indicator at a declared ``(alpha, gamma)`` cell.

    The capacity uses the optimism convention ``L + gamma (U - L)`` and the
    threshold the conservative-opposite convention ``U - gamma (U - L)``.
    """

    capacity_lower, capacity_upper = alpha_cut(a, alpha)
    threshold_lower, threshold_upper = alpha_cut(b, alpha)
    capacity = capacity_lower + gamma * (capacity_upper - capacity_lower)
    threshold = threshold_upper - gamma * (threshold_upper - threshold_lower)
    return 1.0 if capacity >= threshold else 0.0


def proportional_capacity_tfn(modal_capacity: float, left_spread: float,
                              right_spread: float) -> TFN:
    """Return the aggregated capacity TFN under proportional spreads.

    With ``w_i = (w_i (1 - d_l), w_i, w_i (1 + d_r))`` the extension-principle
    sum over any subset is ``(W (1 - d_l), W, W (1 + d_r))`` where ``W`` is the
    modal total.  The success degree therefore depends on the *modal* total
    only, which keeps the dynamic-programming state one-dimensional.  This is
    the property that makes the family computationally free; it does not hold
    for additive or subset-dependent spreads.
    """

    if modal_capacity < 0.0:
        raise ValueError("modal capacity must be non-negative")
    return TFN(modal_capacity * (1.0 - left_spread), modal_capacity,
               modal_capacity * (1.0 + right_spread))


MEASURES: dict[str, Callable[[TFN, TFN], float]] = {
    "pos": possibility_geq,
    "nec": necessity_geq,
    "cr": credibility_geq,
}


def membership_profile(threshold_tfn: TFN, measure: str, left_spread: float,
                       right_spread: float, scale: int, capacity_cap: int,
                       alpha: float | None = None,
                       gamma: float | None = None) -> list[float]:
    """Tabulate the success degree for every scaled modal capacity level.

    ``scale`` is the rational-scaling denominator already used by the solver,
    so entry ``j`` corresponds to modal capacity ``j / scale``.  The returned
    list has length ``capacity_cap + 1`` and is non-decreasing.
    """

    if capacity_cap < 0:
        raise ValueError("capacity cap must be non-negative")
    if scale <= 0:
        raise ValueError("scale must be positive")

    profile: list[float] = []
    for level in range(capacity_cap + 1):
        capacity = proportional_capacity_tfn(level / scale, left_spread, right_spread)
        if measure == "hurwicz":
            if alpha is None or gamma is None:
                raise ValueError("hurwicz measure requires alpha and gamma")
            degree = hurwicz_indicator(capacity, threshold_tfn, alpha, gamma)
        else:
            try:
                degree = MEASURES[measure](capacity, threshold_tfn)
            except KeyError as error:
                raise ValueError(f"unknown measure: {measure}") from error
        profile.append(degree)

    # Monotonicity is a modelling requirement, not an accident: the exact
    # bounds in the manuscript rely on it.  Enforce it defensively so a future
    # membership family cannot silently invalidate those proofs.
    for previous, current in zip(profile, profile[1:]):
        if current < previous - 1e-12:
            raise AssertionError("membership profile must be non-decreasing")
    return profile


def saturation_level(threshold_tfn: TFN, measure: str, left_spread: float,
                     scale: int) -> int:
    """Return the smallest scaled modal capacity whose success degree is 1.

    ``Pos`` saturates at the modal threshold; ``Nec`` saturates when the lower
    capacity endpoint reaches the upper threshold endpoint.  Capping the
    capacity PMF at this level is exact, because no further capacity can raise
    the degree above one.
    """

    if measure == "pos":
        capacity = threshold_tfn.mode
    elif measure in ("nec", "cr"):
        if left_spread >= 1.0:
            raise ValueError("left spread must be below 1 for a finite necessity cap")
        capacity = threshold_tfn.right / (1.0 - left_spread)
    else:
        raise ValueError(f"saturation level is undefined for measure: {measure}")
    level = int(capacity * scale)
    while level / scale < capacity:
        level += 1
    return level
