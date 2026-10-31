"""Independent exhaustive check of selected unit-specific RTS-79 optima (exact rational arithmetic).

Enumerates all 3^10 fleets of the RTS-79 catalog, evaluates each with unit-specific Hurwicz
capacities and the exact reliability (all up/down states of the installed units), and applies the
lexicographic rule (minimum cost, then maximum reliability).  It shares no code with the C++
search used by run_power_rts_unit_specific_v7.py, and checks four results of that experiment:

  * the common-spread optimum at (0.3, 0.7)              -> 158
  * the unit-specific optima at (0, 0), (0.3, 0), (0.3, 0.7) -> 203, 180, 180
  * the 171 M$ fleet (1 x 155 MW coal + 2 x 20 MW CT) selected at (0.3, 0.7) when the CT
    downside is 0.10, and at (0.5, 0.7) in the hydro-only variant, and that this fleet is the
    crisp optimum of no integer threshold of the parametric solution.

Writes results/generated/unit_specific_bruteforce_v8.json.

    PYTHONPATH=src python3 experiments/verify_unit_specific_bruteforce_v8.py
"""

from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS  # noqa: E402

B = BENCHMARKS["power_rts_adapted"]
RESULTS = ROOT / "results" / "generated"

COMMON = {mw: (0.10, 0.10) for mw in (12, 20, 50, 76, 100, 155)}
UNIT = {12: (0.05, 0.02), 20: (0.15, 0.05), 50: (0.10, 0.0), 76: (0.05, 0.02), 100: (0.05, 0.02),
        155: (0.05, 0.02)}
CT_010 = {**UNIT, 20: (0.10, 0.05)}
HYDRO_ONLY = {**COMMON, 50: (0.10, 0.0)}

CHECKS = [  # (label, alpha, gamma, spreads, expected cost)
    ("common spread", 0.3, 0.7, COMMON, 158),
    ("unit-specific", 0.0, 0.0, UNIT, 203),
    ("unit-specific", 0.3, 0.0, UNIT, 180),
    ("unit-specific", 0.3, 0.7, UNIT, 180),
    ("unit-specific, CT delta_l = 0.10", 0.3, 0.7, CT_010, 171),
    ("hydro only", 0.5, 0.7, HYDRO_ONLY, 171),
]


def solve(alpha: float, gamma: float, spreads: dict) -> tuple[int, F, tuple[int, ...]]:
    a, g = F(str(alpha)), F(str(gamma))
    k1, k2, k3 = (F(str(x)) for x in B["K_tfn"])
    lower, upper = k1 + a * (k2 - k1), k3 - a * (k3 - k2)
    threshold = upper - g * (upper - lower)

    def capacity(mw: int) -> F:
        left, right = (F(str(x)) for x in spreads[mw])
        w = F(mw)
        lo = w * (1 - left) + a * left * w
        hi = w * (1 + right) - a * right * w
        return lo + g * (hi - lo)

    target = F(str(B["R_target"]))
    best = None
    for design in itertools.product(*[range(len(c)) for c in B["choices"]]):
        units = [B["choices"][i][c] for i, c in enumerate(design) if c]
        cost = sum(u[1] for u in units)
        if best is not None and cost > best[0]:
            continue
        caps = [capacity(u[2]) for u in units]
        if sum(caps) < threshold:
            continue
        rel = F(0)
        for states in itertools.product((0, 1), repeat=len(units)):
            if sum(c for c, s in zip(caps, states) if s) >= threshold:
                p = F(1)
                for u, s in zip(units, states):
                    pu = F(str(u[0]))
                    p *= pu if s else 1 - pu
                rel += p
        if rel < target:
            continue
        if best is None or cost < best[0] or (cost == best[0] and rel > best[1]):
            best = (cost, rel, design)
    return best


def main() -> None:
    library = json.loads((RESULTS / "parametric_stability_rts.json").read_text())["parametric_path"]
    crisp = {tuple(e["chromosome"]) for e in library if e.get("chromosome")}
    out = []
    for label, a, g, spreads, expected in CHECKS:
        cost, rel, design = solve(a, g, spreads)
        row = {"check": label, "alpha": a, "gamma": g, "cost": cost, "reliability": float(rel),
               "design": list(design), "expected_cost": expected, "agrees": cost == expected,
               "in_crisp_library": tuple(design) in crisp}
        out.append(row)
        print(row)
        assert cost == expected, row
    outside = [r for r in out if r["cost"] == 171]
    assert outside and all(not r["in_crisp_library"] for r in outside)
    (RESULTS / "unit_specific_bruteforce_v8.json").write_text(json.dumps(
        {"description": __doc__.split("\n\n")[0], "checks": out}, indent=2) + "\n")


if __name__ == "__main__":
    main()
