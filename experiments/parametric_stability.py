"""Exact parametric solution and stability regions over the whole (alpha, gamma) square.

Dimensional reduction
---------------------
Under a common proportional spread ``d`` every scalarised capacity is
``w_i * m`` with ``m = 1 + d (1 - alpha) (2 gamma - 1)``.  The success event
``sum_i w_i m x_i >= k_H`` is therefore equivalent to

    sum_i w_i x_i >= kappa(alpha, gamma),   kappa = k_H / m,

so the two-parameter family acts on the design problem through the single
scalar ``kappa`` -- an *effective modal threshold*.  Both adapted catalogs have
integer modal capacities, hence the success event depends only on the integer
``T = ceil(kappa)``.

Consequences used here:

* the complete parametric solution is obtained by solving one crisp problem per
  attainable ``T``, which covers the whole continuum ``[0,1]^2`` rather than a
  finite grid sample;
* the optimal design is piecewise constant in ``(alpha, gamma)``, and the
  region boundaries are the closed-form level sets ``kappa(alpha, gamma) = T``.

Run from the repository root:

    PYTHONPATH=src python3 experiments/parametric_stability.py
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import ceil
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.solver import solve_exact

RESULTS = ROOT / "results" / "generated"


def kappa(threshold_tfn: tuple[float, float, float], alpha: Fraction, gamma: Fraction,
          spread: Fraction) -> Fraction:
    """Effective modal threshold at a declared cell (exact rational arithmetic)."""

    k1, k2, k3 = (Fraction(str(value)) for value in threshold_tfn)
    cut_lower = k1 + alpha * (k2 - k1)
    cut_upper = k3 - alpha * (k3 - k2)
    demand = cut_upper - gamma * (cut_upper - cut_lower)
    factor = 1 + spread * (1 - alpha) * (2 * gamma - 1)
    if factor <= 0:
        raise ValueError("capacity factor must stay positive")
    return demand / factor


def kappa_range(threshold_tfn, spread: Fraction) -> tuple[Fraction, Fraction]:
    """kappa is largest at (0,0) and smallest at (0,1); verified by enumeration."""

    corners = [kappa(threshold_tfn, Fraction(a), Fraction(g), spread)
               for a in (0, 1) for g in (0, 1)]
    return min(corners), max(corners)


def solve_at_integer_threshold(benchmark: dict, threshold: int, backend: str = "python") -> dict:
    """Exact lexicographic solution of the crisp problem with modal capacities.

    The C++ backend is a faithful port of the Python search and is validated
    cell-by-cell against the archived audits before use; it changes runtime,
    not the certified optimum.
    """

    if backend == "cpp":
        from fuzzy_wkon.cpp_backend import solve_at_integer_threshold_cpp
        return solve_at_integer_threshold_cpp(benchmark, threshold)
    crisp = dict(benchmark)
    crisp["K_tfn"] = (float(threshold), float(threshold), float(threshold))
    return solve_exact(crisp, alpha=0.0, gamma=0.0, left_spread=0.0, right_spread=0.0)


def run(spread_value: float = 0.1, only: str | None = None, verbose: bool = False,
        backend: str = "python") -> dict:
    spread = Fraction(str(spread_value))
    report: dict = {"spread": spread_value, "backend": backend, "catalogs": {}}

    for key, benchmark in BENCHMARKS.items():
        if only is not None and key != only:
            continue
        low, high = kappa_range(benchmark["K_tfn"], spread)
        thresholds = list(range(ceil(low), ceil(high) + 1))

        path = []
        for threshold in thresholds:
            started = perf_counter()
            solution = solve_at_integer_threshold(benchmark, threshold, backend)
            if verbose:
                print(f"    [{key}] T={threshold} solved in {perf_counter()-started:.1f}s "
                      f"cost={solution.get('best_cost')}", flush=True)
            path.append({
                "T": threshold,
                "feasible": bool(solution.get("feasible")),
                "cost": solution.get("best_cost"),
                "reliability": solution.get("best_reliability"),
                "chromosome": solution.get("best_chromosome"),
            })

        # Stability regions: maximal runs of T sharing one optimal design.
        regions = []
        for entry in path:
            signature = ("infeasible" if not entry["feasible"]
                         else "-".join(map(str, entry["chromosome"])))
            if regions and regions[-1]["signature"] == signature:
                regions[-1]["T_max"] = entry["T"]
            else:
                regions.append({
                    "signature": signature,
                    "T_min": entry["T"],
                    "T_max": entry["T"],
                    "feasible": entry["feasible"],
                    "cost": entry["cost"],
                    "reliability": entry["reliability"],
                    "chromosome": entry["chromosome"],
                })

        report["catalogs"][key] = {
            "kappa_min": float(low),
            "kappa_max": float(high),
            "integer_thresholds": thresholds,
            "parametric_path": path,
            "stability_regions": regions,
            "n_distinct_designs": len({r["signature"] for r in regions if r["feasible"]}),
            "n_breakpoints": max(0, len(regions) - 1),
        }

    return report


def verify_against_audit(report: dict, spread_value: float = 0.1) -> dict:
    """Every archived V2 audit cell must be recovered by the parametric path."""

    spread = Fraction(str(spread_value))
    audits = {
        "pharma_adapted": RESULTS / "pharma_oracle_audit.json",
        "fyffe_adapted": RESULTS / "fyffe_exact_audit.json",
        "power_rts_adapted": RESULTS / "benchmark_power_rts.json",
    }
    summary = {}
    for key, benchmark in BENCHMARKS.items():
        raw = json.loads(audits[key].read_text())
        cells = raw if isinstance(raw, list) else raw["results"]
        if key not in report["catalogs"]:
            continue
        lookup = {entry["T"]: entry for entry in report["catalogs"][key]["parametric_path"]}
        checked = matched_cost = matched_rel = matched_design = 0
        mismatches = []
        for cell in cells:
            solver = cell.get("solver", cell)
            archived_cost = solver.get("best_cost", solver.get("cost"))
            archived_rel = solver.get("best_reliability", solver.get("best_reliability_native", solver.get("R_hurwicz")))
            archived_chrom = solver.get("best_chromosome", solver.get("chromosome", solver.get("design")))
            threshold = ceil(kappa(benchmark["K_tfn"], Fraction(str(cell["alpha"])),
                                   Fraction(str(cell["gamma"])), spread))
            predicted = lookup.get(threshold)
            if predicted is None:
                continue
            checked += 1
            cost_ok = (predicted["cost"] == archived_cost)
            rel_ok = (not predicted["feasible"] or abs(predicted["reliability"] - archived_rel) < 1e-4)
            design_ok = (not predicted["feasible"] or predicted["chromosome"] == archived_chrom)

            if cost_ok:
                matched_cost += 1
            if rel_ok:
                matched_rel += 1
            if design_ok:
                matched_design += 1

            if not (cost_ok and rel_ok and design_ok):
                mismatches.append({"alpha": cell["alpha"], "gamma": cell["gamma"],
                                   "T": threshold, "archived_cost": archived_cost,
                                   "parametric_cost": predicted["cost"],
                                   "cost_ok": cost_ok, "rel_ok": rel_ok, "design_ok": design_ok})
        summary[key] = {
            "checked": checked,
            "matched": matched_cost,
            "matched_cost": matched_cost,
            "matched_reliability": matched_rel,
            "matched_design": matched_design,
            "all_matched": (matched_cost == checked and matched_rel == checked and matched_design == checked),
            "mismatches": mismatches,
        }
    return summary


if __name__ == "__main__":
    import sys

    arguments = [a for a in sys.argv[1:] if not a.startswith("--")]
    backend = "cpp" if "--cpp" in sys.argv[1:] else "python"
    only = arguments[0] if arguments else None
    report = run(only=only, verbose=True, backend=backend)
    report["audit_recovery"] = verify_against_audit(report)
    suffix = f"_{only}" if only else ""
    output = RESULTS / f"parametric_stability{suffix}.json"
    output.write_text(json.dumps(report, indent=2))
    if "power_rts_adapted" in report["catalogs"]:
        rts_output = RESULTS / "parametric_stability_rts.json"
        rts_output.write_text(json.dumps(report["catalogs"]["power_rts_adapted"], indent=2))
        print(f"written: {rts_output}")

    for key, block in report["catalogs"].items():
        print(f"\n[{key}]  kappa in [{block['kappa_min']:.3f}, {block['kappa_max']:.3f}]"
              f"  -> {len(block['integer_thresholds'])} distinct crisp problems")
        print(f"  distinct optimal designs: {block['n_distinct_designs']}"
              f"   breakpoints: {block['n_breakpoints']}")
        for region in block["stability_regions"]:
            label = "INFEASIBLE" if not region["feasible"] else f"cost={region['cost']:>4}"
            print(f"    T in [{region['T_min']:>3}, {region['T_max']:>3}]  {label}"
                  f"   {region['signature']}")
        recovery = report["audit_recovery"][key]
        print(f"  audit recovery: {recovery['matched_cost']}/{recovery['checked']} cost, "
              f"{recovery['matched_reliability']}/{recovery['checked']} rel, "
              f"{recovery['matched_design']}/{recovery['checked']} design reproduced"
              f"  (mismatches: {len(recovery['mismatches'])})")
    print(f"\nwritten: {output}")
