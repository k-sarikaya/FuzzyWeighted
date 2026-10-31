"""Explicit fuzzy interval and Hurwicz decision semantics for V2.

The module deliberately does not call the resulting reliability
"possibilistic".  It implements a named interval decision rule and leaves
possibility/necessity measures as separate alternatives.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class TFN:
    left: float
    mode: float
    right: float

    def __post_init__(self) -> None:
        if not (self.left <= self.mode <= self.right):
            raise ValueError("TFN must satisfy left <= mode <= right")

    def cut(self, alpha: float) -> Tuple[float, float]:
        return alpha_cut(self, alpha)


def _check_unit(value: float, name: str) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must belong to [0, 1]")


def alpha_cut(tfn: TFN, alpha: float) -> Tuple[float, float]:
    _check_unit(alpha, "alpha")
    lower = tfn.left + alpha * (tfn.mode - tfn.left)
    upper = tfn.right - alpha * (tfn.right - tfn.mode)
    return lower, upper


def hurwicz_value(tfn: TFN, alpha: float, gamma: float) -> float:
    """Return the Hurwicz scalar for an alpha-cut interval.

    gamma=0 selects the lower endpoint and gamma=1 selects the upper
    endpoint.  For a threshold the caller can use ``1-gamma`` when a larger
    threshold is more conservative.
    """

    _check_unit(gamma, "gamma")
    lower, upper = alpha_cut(tfn, alpha)
    return lower + gamma * (upper - lower)


def threshold_hurwicz_value(tfn: TFN, alpha: float, gamma: float) -> float:
    """Return the demand threshold under the model's optimism convention."""

    _check_unit(gamma, "gamma")
    lower, upper = alpha_cut(tfn, alpha)
    return upper - gamma * (upper - lower)


def alpha_derivative(tfn: TFN, gamma: float) -> float:
    """Derivative of the Hurwicz capacity scalar with respect to alpha."""

    _check_unit(gamma, "gamma")
    return (1.0 - gamma) * (tfn.mode - tfn.left) - gamma * (tfn.right - tfn.mode)


def gamma_from_certainty_equivalent(lower: float, upper: float, certainty_equivalent: float) -> float:
    """Map an elicited interval certainty equivalent to gamma."""

    if upper <= lower:
        raise ValueError("upper must be greater than lower")
    if not lower <= certainty_equivalent <= upper:
        raise ValueError("certainty equivalent must lie within [lower, upper]")
    return (certainty_equivalent - lower) / (upper - lower)
