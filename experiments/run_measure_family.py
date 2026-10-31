"""Evaluate the possibility/necessity/Hurwicz reliability family on both catalogs.

Two things are established here.

1. *Regression gate.*  Every Hurwicz reliability recomputed through the graded
   kernel must reproduce the archived V2 audit value.  The decision-rule family
   is an extension, never a revision, of the verified V2 numbers.

2. *Decomposition identity.*  For each selected design the expected-necessity
   and expected-possibility reliabilities are compared with the alpha-integral
   of the ``gamma=0`` and ``gamma=1`` columns of the Hurwicz grid.

Run from the repository root:

    PYTHONPATH=src python3 experiments/run_measure_family.py
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import lcm
from pathlib import Path

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.fuzzy import TFN
from fuzzy_wkon.membership import membership_profile, saturation_level
from fuzzy_wkon.reliability import measure_reliability
from fuzzy_wkon.solver import prepare_options

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "generated"
AUDITS = {
    "pharma_adapted": RESULTS / "pharma_oracle_audit.json",
    "fyffe_adapted": RESULTS / "fyffe_exact_audit.json",
}
ALPHA_QUADRATURE = 20000


def modal_options(benchmark: dict, chromosome: list[int]) -> tuple[list[tuple[float, int]], int]:
    """Return selected options with *modal* integer capacities and their scale.

    The graded family is indexed by the modal total capacity, so the Hurwicz
    scalarisation must not be baked into the weights here.
    """

    denominators = [Fraction(str(value)).denominator for value in benchmark["K_tfn"]]
    for choices in benchmark["choices"]:
        denominators.extend(Fraction(str(weight)).denominator for _, _, weight in choices)
    scale = 1
    for denominator in denominators:
        scale = lcm(scale, denominator)

    selected: list[tuple[float, int]] = []
    for choices, choice_index in zip(benchmark["choices"], chromosome):
        reliability, _cost, weight = choices[choice_index]
        if choice_index == 0:
            continue
        selected.append((float(reliability), int(Fraction(str(weight)) * scale)))
    return selected, scale


def hurwicz_profile(threshold: TFN, alpha: float, gamma: float, spread: float,
                    scale: int, cap: int) -> list[float]:
    profile: list[float] = []
    for level in range(cap + 1):
        modal = level / scale
        lower = modal * (1.0 - spread) + alpha * spread * modal
        upper = modal * (1.0 + spread) - alpha * spread * modal
        capacity = lower + gamma * (upper - lower)
        cut_lower = threshold.left + alpha * (threshold.mode - threshold.left)
        cut_upper = threshold.right - alpha * (threshold.right - threshold.mode)
        demand = cut_upper - gamma * (cut_upper - cut_lower)
        profile.append(1.0 if capacity >= demand else 0.0)
    profile[-1] = 1.0
    return profile


def alpha_integral(options, threshold: TFN, gamma: float, spread: float,
                   scale: int, cap: int, points: int = ALPHA_QUADRATURE) -> float:
    total = 0.0
    for index in range(points):
        alpha = (index + 0.5) / points
        total += measure_reliability(options, hurwicz_profile(threshold, alpha, gamma, spread, scale, cap))
    return total / points


def run(spread: float = 0.1) -> dict:
    report: dict = {"spread": spread, "alpha_quadrature": ALPHA_QUADRATURE, "catalogs": {}}

    for key, benchmark in BENCHMARKS.items():
        if key not in AUDITS:
            continue
        audit = json.loads(AUDITS[key].read_text())
        cells = audit if isinstance(audit, list) else audit["results"]
        threshold = TFN(*benchmark["K_tfn"])
        rows = []
        worst_regression = 0.0

        for cell in cells:
            solver = cell.get("solver", cell)
            chromosome = solver.get("best_chromosome") or solver.get("chromosome")
            if chromosome is None or solver.get("feasible") is False:
                rows.append({"alpha": cell["alpha"], "gamma": cell["gamma"], "feasible": False})
                continue
            options, scale = modal_options(benchmark, chromosome)
            cap = max(saturation_level(threshold, "pos", spread, scale),
                      saturation_level(threshold, "nec", spread, scale))

            recomputed = measure_reliability(
                options,
                hurwicz_profile(threshold, cell["alpha"], cell["gamma"], spread, scale, cap),
            )
            archived = solver.get("best_reliability")
            if archived is None:
                archived = solver["best_reliability_python"]
            worst_regression = max(worst_regression, abs(recomputed - archived))

            rows.append({
                "alpha": cell["alpha"],
                "gamma": cell["gamma"],
                "feasible": True,
                "cost": solver["best_cost"],
                "chromosome": chromosome,
                "R_hurwicz_archived": archived,
                "R_hurwicz_recomputed": recomputed,
                "R_possibility": measure_reliability(
                    options, membership_profile(threshold, "pos", spread, spread, scale, cap)),
                "R_necessity": measure_reliability(
                    options, membership_profile(threshold, "nec", spread, spread, scale, cap)),
            })

        # Decomposition identity on the three distinct designs with the widest
        # cost spread, so the check is not dominated by one repeated chromosome.
        feasible = [row for row in rows if row["feasible"]]
        unique: dict[str, dict] = {}
        for row in feasible:
            unique.setdefault("-".join(map(str, row["chromosome"])), row)
        probes = sorted(unique.values(), key=lambda row: row["cost"])
        probes = [probes[0], probes[len(probes) // 2], probes[-1]] if len(probes) >= 3 else probes

        identity = []
        for row in probes:
            options, scale = modal_options(benchmark, row["chromosome"])
            cap = max(saturation_level(threshold, "pos", spread, scale),
                      saturation_level(threshold, "nec", spread, scale))
            identity.append({
                "chromosome": row["chromosome"],
                "cost": row["cost"],
                "R_necessity": row["R_necessity"],
                "integral_gamma_0": alpha_integral(options, threshold, 0.0, spread, scale, cap),
                "R_possibility": row["R_possibility"],
                "integral_gamma_1": alpha_integral(options, threshold, 1.0, spread, scale, cap),
            })
            identity[-1]["necessity_gap"] = abs(identity[-1]["R_necessity"] - identity[-1]["integral_gamma_0"])
            identity[-1]["possibility_gap"] = abs(identity[-1]["R_possibility"] - identity[-1]["integral_gamma_1"])

        report["catalogs"][key] = {
            "cells": rows,
            "worst_hurwicz_regression": worst_regression,
            "decomposition_identity": identity,
        }

    return report


if __name__ == "__main__":
    report = run()
    output = RESULTS / "measure_family.json"
    output.write_text(json.dumps(report, indent=2))
    for key, block in report["catalogs"].items():
        print(f"[{key}] worst Hurwicz regression vs archived audit: {block['worst_hurwicz_regression']:.3e}")
        for row in block["decomposition_identity"]:
            print(f"    cost={row['cost']:>4}  |R_nec - int(gamma=0)|={row['necessity_gap']:.3e}"
                  f"   |R_pos - int(gamma=1)|={row['possibility_gap']:.3e}")
    print(f"written: {output}")
