"""Python driver for the C++ exact backend (`cpp/exact_solver.cpp`).

The audited Python layer keeps ownership of the exact rational scalarisation;
the C++ program only performs the search on the resulting *integer* problem.
Two consequences matter for the manuscript:

* no fuzzy arithmetic is duplicated across languages, so the two backends
  cannot silently disagree about the model;
* the gcd reduction applied here is the value-preserving rational scaling of
  Proposition 1, so it changes runtime but not the certified optimum.

`solve_exact_cpp` mirrors the return shape of `solver.solve_exact` and adds a
``status`` field so that a censored (timed-out) observation is never reported
as an exact optimum.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import os
from functools import lru_cache
from math import gcd
from pathlib import Path
from typing import Any

from .solver import prepare_options

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "cpp" / "exact_solver.cpp"
BINARY = ROOT / "cpp" / ("exact_solver.exe" if os.name == "nt" else "exact_solver")


class BackendUnavailable(RuntimeError):
    """Raised when the C++ backend can neither be found nor built."""


@lru_cache(maxsize=1)
def ensure_binary() -> Path:
    """Build the backend once per process if it is missing or out of date."""

    if BINARY.exists() and BINARY.stat().st_mtime >= SOURCE.stat().st_mtime:
        return BINARY
    compiler = shutil.which("g++") or shutil.which("clang++")
    if compiler is None:
        raise BackendUnavailable("no C++ compiler (g++ or clang++) on PATH")
    completed = subprocess.run(
        [compiler, "-O3", "-std=c++17", str(SOURCE), "-o", str(BINARY)],
        capture_output=True, text=True,
    )
    if completed.returncode != 0:
        raise BackendUnavailable(f"compilation failed:\n{completed.stderr}")
    return BINARY


def build_payload(options, threshold: int, target: float, time_limit: float) -> tuple[str, int]:
    """Serialise an already-integral problem, reduced by the common divisor."""

    divisor = threshold
    for row in options:
        for option in row:
            divisor = gcd(divisor, option.weight_scaled)
    divisor = divisor or 1

    lines = [f"{len(options)} {target!r} {threshold // divisor} {time_limit!r}"]
    for row in options:
        lines.append(str(len(row)))
        for option in row:
            lines.append(f"{option.choice_index} {option.reliability!r} "
                         f"{option.cost} {option.weight_scaled // divisor}")
    return "\n".join(lines) + "\n", divisor


def solve_integer_problem(options, threshold: int, target: float,
                          time_limit: float = 0.0) -> dict[str, Any]:
    binary = ensure_binary()
    payload, divisor = build_payload(options, threshold, target, time_limit)
    completed = subprocess.run([str(binary)], input=payload, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(f"backend failed: {completed.stderr.strip()}")
    result = json.loads(completed.stdout.strip())
    result["common_divisor"] = divisor
    return result


def solve_exact_cpp(benchmark: dict[str, Any], alpha: float, gamma: float,
                    left_spread: float = 0.1, right_spread: float = 0.1,
                    time_limit: float = 0.0,
                    unit_spreads: list[list[tuple[float, float]]] | None = None) -> dict[str, Any]:
    """Exact lexicographic solution through the C++ backend."""

    options, threshold, scale = prepare_options(benchmark, alpha, gamma,
                                                left_spread, right_spread,
                                                unit_spreads=unit_spreads)
    result = solve_integer_problem(options, threshold, benchmark["R_target"], time_limit)
    result.update({"alpha": alpha, "gamma": gamma, "scale": scale,
                   "threshold_scaled": threshold, "backend": "cpp"})
    return result


def solve_at_integer_threshold_cpp(benchmark: dict[str, Any], threshold: int,
                                   time_limit: float = 0.0) -> dict[str, Any]:
    """Crisp modal problem at an effective integer threshold (Theorem 1 of the manuscript)."""

    crisp = dict(benchmark)
    crisp["K_tfn"] = (float(threshold), float(threshold), float(threshold))
    result = solve_exact_cpp(crisp, alpha=0.0, gamma=0.0,
                             left_spread=0.0, right_spread=0.0, time_limit=time_limit)
    result["T"] = threshold
    return result
