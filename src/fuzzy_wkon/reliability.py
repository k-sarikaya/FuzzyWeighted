"""Exact reliability kernels used by both the solver and independent tests."""

from __future__ import annotations

from itertools import product
from typing import Iterable, Mapping, Sequence


def capped_distribution(options: Iterable[tuple[float, int]], threshold: int) -> list[float]:
    """Convolve independent options and cap all weight >= threshold."""

    if threshold < 0:
        raise ValueError("threshold must be non-negative")
    pmf = [0.0] * (threshold + 1)
    pmf[0] = 1.0
    for probability, weight in options:
        if weight < 0:
            raise ValueError("weights must be non-negative")
        nxt = [0.0] * (threshold + 1)
        for current, mass in enumerate(pmf):
            if not mass:
                continue
            nxt[current] += mass * (1.0 - probability)
            nxt[min(threshold, current + weight)] += mass * probability
        pmf = nxt
    return pmf


def static_reliability(options: Iterable[tuple[float, int]], threshold: int) -> float:
    return capped_distribution(options, threshold)[threshold]


def subset_enumeration_reliability(options: Sequence[tuple[float, int]], threshold: int) -> float:
    """Independent oracle; intentionally exponential and used only in tests."""

    result = 0.0
    for states in product((0, 1), repeat=len(options)):
        weight = sum(option[1] for option, state in zip(options, states) if state)
        if weight < threshold:
            continue
        probability = 1.0
        for state, (p, _) in zip(states, options):
            probability *= p if state else (1.0 - p)
        result += probability
    return result


def distribution_reliability(pmf: Sequence[float]) -> float:
    return float(pmf[-1])


def measure_reliability(options: Iterable[tuple[float, int]],
                        membership_profile: Sequence[float]) -> float:
    """Probability of the fuzzy success event under a graded membership.

    ``membership_profile[j]`` is the success degree at scaled modal capacity
    ``j``; the last entry must equal one so that capping the convolution is
    exact.  With an indicator profile this reduces to ``static_reliability``,
    which is why every V2 number remains reproducible inside the family.
    """

    cap = len(membership_profile) - 1
    if cap < 0:
        raise ValueError("membership profile must be non-empty")
    if abs(membership_profile[-1] - 1.0) > 1e-12:
        raise ValueError("membership profile must saturate at one")
    pmf = capped_distribution(options, cap)
    return float(sum(mass * degree for mass, degree in zip(pmf, membership_profile)))


def subset_enumeration_measure_reliability(options: Sequence[tuple[float, int]],
                                           membership_profile: Sequence[float]) -> float:
    """Independent exponential oracle for :func:`measure_reliability`."""

    cap = len(membership_profile) - 1
    result = 0.0
    for states in product((0, 1), repeat=len(options)):
        weight = sum(option[1] for option, state in zip(options, states) if state)
        degree = membership_profile[min(cap, weight)]
        if not degree:
            continue
        probability = 1.0
        for state, (p, _) in zip(states, options):
            probability *= p if state else (1.0 - p)
        result += probability * degree
    return result
