"""Deterministic synthetic catalog families for bounded scaling evidence."""

from __future__ import annotations

import random


PROFILES = {
    "screening": {"threshold_ratio": 0.34, "reliability_target": 0.80},
    "design": {"threshold_ratio": 0.48, "reliability_target": 0.90},
    "stringent": {"threshold_ratio": 0.60, "reliability_target": 0.95},
}


def synthetic_benchmark(n: int, profile: str, seed: int = 20260812) -> dict:
    """Create a fully declared three-choice catalog (null/base/premium).

    The generator is evidence for computational profiling only. It is not an
    industrial data model and it does not reuse either adapted benchmark.
    """
    if n < 2:
        raise ValueError("n must be at least two")
    if profile not in PROFILES:
        raise ValueError(f"unknown profile: {profile}")
    rng = random.Random(seed + 1009 * n + 7919 * list(PROFILES).index(profile))
    choices = []
    maximum_capacity = 0
    for _ in range(n):
        base_weight = rng.randint(5, 9)
        premium_weight = base_weight + rng.randint(2, 5)
        base_reliability = round(rng.uniform(0.78, 0.90), 4)
        premium_reliability = round(min(0.985, base_reliability + rng.uniform(0.04, 0.09)), 4)
        base_cost = rng.randint(5, 10)
        premium_cost = base_cost + rng.randint(3, 7)
        choices.append([
            (0.0, 0, 0),
            (base_reliability, base_cost, base_weight),
            (premium_reliability, premium_cost, premium_weight),
        ])
        maximum_capacity += premium_weight
    settings = PROFILES[profile]
    modal_threshold = max(1, round(settings["threshold_ratio"] * maximum_capacity))
    return {
        "name": f"deterministic scaling family n={n} profile={profile}",
        "N": n,
        # Round the declared decimal endpoints before Fraction(str(value)) is
        # applied by the exact solver. This prevents binary floating-point
        # artifacts (for example 56.10000000000001) from creating a spurious
        # common denominator and an artificial state-space explosion.
        "K_tfn": (round(0.9 * modal_threshold, 6), float(modal_threshold), round(1.1 * modal_threshold, 6)),
        "R_target": settings["reliability_target"],
        "choices": choices,
        "generator": {
            "seed": seed,
            "profile": profile,
            **settings,
            "options_per_subsystem": 3,
        },
    }
