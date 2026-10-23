"""Unit-specific capacity spreads on the IEEE RTS-79 generation fleet.

The baseline RTS-79 study uses one proportional spread (delta = 0.10) for every
unit, which is what Theorem 1 requires.  This experiment replaces it with
unit-specific spreads and re-solves every policy cell with the unreduced exact
search, because with heterogeneous spreads the effective threshold is no longer
a closed-form function of (alpha, gamma).

The hydro spread is taken from RTS-79 itself: Table 6 of the 1979 paper gives the
available hydro capacity as 100% in quarters 1-2 and 90% in quarters 3-4, so the
hydro capacity TFN is (45, 50, 50) MW, i.e. (delta_l, delta_r) = (0.10, 0).  The
thermal spreads are an assumption of this study (RTS-79 Table 4 covers full outages
only and leaves partial outages to the user): combustion turbines strongly
sensitive to ambient temperature, steam units mildly sensitive.  The hydro sweep
below varies the downside delta_l of hydro to show what the Table 6 value implies
relative to a temperature-only view (delta_l = 0).

What is recorded, per policy cell:
  * common-spread and unit-specific optima, cost, fleet, unreduced runtime;
  * whether the unit-specific optimum belongs to the crisp library (the
    optima of the 73 integer thresholds of the parametric solution) and at
    which thresholds;
  * for every cell whose optimum changes, a 2x2 decomposition
    (common / thermal-only / hydro-only / unit-specific) and two
    sensitivity sweeps (delta_hydro; CT and steam spreads).

Produces:
  results/generated/rts_unit_specific_v7.json
  results/generated/table_rts_unit_specific.tex, manuscript/table_rts_unit_specific.tex

Run from the repository root (about 1-1.5 h on two cores):
    PYTHONPATH=src python3 experiments/run_power_rts_unit_specific_v7.py
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS                      # noqa: E402
from fuzzy_wkon.cpp_backend import ensure_binary, solve_exact_cpp  # noqa: E402

ALPHAS = (0.0, 0.3, 0.5, 0.7, 0.9, 1.0)
GAMMAS = (0.0, 0.3, 0.5, 0.7, 1.0)
COMMON = 0.10
TIME_LIMIT = 600.0
HYDRO_MW = 50
HYDRO_SWEEP = (0.0, 0.02, 0.05, 0.10)   # hydro delta_l; delta_r stays 0 (Table 6: no capacity above 100%)
THERMAL_PERTURBATIONS = (
    ("base", 0.15, 0.05, 0.05, 0.02),
    ("CT delta_l 0.20", 0.20, 0.05, 0.05, 0.02),
    ("CT delta_l 0.10", 0.10, 0.05, 0.05, 0.02),
    ("steam delta_l 0.08", 0.15, 0.05, 0.08, 0.02),
    ("steam delta_l 0.02", 0.15, 0.05, 0.02, 0.02),
)
BENCH = "power_rts_adapted"
NAMES = {(12, 14): "12 MW oil", (20, 18): "20 MW CT", (50, 45): "50 MW hydro",
         (76, 68): "76 MW coal", (100, 85): "100 MW oil", (155, 135): "155 MW coal"}
ORDER = ["50 MW hydro", "76 MW coal", "155 MW coal", "100 MW oil", "12 MW oil", "20 MW CT"]


def load_unit_spreads() -> dict[int, tuple[float, float]]:
    with open(ROOT / "data" / "original" / "ieee_rts79_units.csv", encoding="utf-8") as f:
        return {int(r["capacity_mw"]): (float(r["derate_low"]), float(r["derate_high"]))
                for r in csv.DictReader(f)}


def spread_structure(benchmark: dict, spreads: dict[int, tuple[float, float]]):
    rows = []
    for bay in benchmark["choices"]:
        row = [(0.0, 0.0)]
        for _p, _c, w in bay[1:]:
            row.append(spreads[w])
        rows.append(row)
    return rows


def describe(benchmark: dict, chromosome: list[int] | None) -> str:
    if not chromosome:
        return "none"
    counts = Counter()
    for bay, c in enumerate(chromosome):
        if c:
            _, cost, mw = benchmark["choices"][bay][c]
            counts[NAMES.get((mw, cost), f"{mw} MW")] += 1
    parts = [f"{counts[u]} $\\times$ {u}" for u in ORDER if u in counts]
    parts += [f"{n} $\\times$ {u}" for u, n in sorted(counts.items()) if u not in ORDER]
    return " + ".join(parts)


def variant_spreads(name: str, unit: dict[int, tuple[float, float]]) -> dict[int, tuple[float, float]]:
    if name == "common":
        return {mw: (COMMON, COMMON) for mw in unit}
    if name == "thermal_only":          # thermal units unit-specific, hydro at the common spread
        s = dict(unit); s[HYDRO_MW] = (COMMON, COMMON); return s
    if name == "hydro_only":            # hydro unit-specific, thermal units at the common spread
        s = {mw: (COMMON, COMMON) for mw in unit}; s[HYDRO_MW] = unit[HYDRO_MW]; return s
    if name == "unit_specific":
        return dict(unit)
    raise ValueError(name)


def solve(job: tuple) -> dict:
    """job = (tag, alpha, gamma, spreads or None). None = common spread, reduced path."""
    tag, alpha, gamma, spreads = job
    benchmark = BENCHMARKS[BENCH]
    t0 = time.time()
    if spreads is None:
        res = solve_exact_cpp(benchmark, alpha, gamma, left_spread=COMMON, right_spread=COMMON)
    else:
        res = solve_exact_cpp(benchmark, alpha, gamma, time_limit=TIME_LIMIT,
                              unit_spreads=spread_structure(benchmark, spreads))
    elapsed = time.time() - t0
    done = res.get("status", "completed") == "completed"
    chrom = res.get("best_chromosome") if done else None
    out = {"tag": tag, "alpha": alpha, "gamma": gamma,
           "status": "completed" if done else res.get("status", "timeout"),
           "cost": res.get("best_cost") if done else None,
           "reliability": res.get("best_reliability") if done else None,
           "chromosome": chrom, "fleet": describe(benchmark, chrom),
           "nodes": res.get("nodes"), "runtime_seconds": elapsed}
    print(f"  {tag} ({alpha},{gamma}): {out['status']} cost={out['cost']} {elapsed:.1f}s", flush=True)
    return out


def crisp_library() -> tuple[dict[tuple, list[int]], int]:
    data = json.loads((ROOT / "results" / "generated" / "parametric_stability_rts.json").read_text())
    library: dict[tuple, list[int]] = {}
    n = 0
    for entry in data.get("parametric_path", []):
        if entry.get("chromosome") is not None:
            n += 1
            library.setdefault(tuple(entry["chromosome"]), []).append(entry["T"])
    return library, n


def run(workers: int) -> dict:
    unit = load_unit_spreads()
    library, n_thresholds = crisp_library()
    cells = [(a, g) for g in GAMMAS for a in ALPHAS]
    ensure_binary()  # build once before forking, so workers do not race on the compiler output

    with Pool(workers) as pool:
        print("Phase 1: common and unit-specific optima in all 30 cells", flush=True)
        jobs = [("common", a, g, None) for a, g in cells] + \
               [("unit_specific", a, g, unit) for a, g in cells]
        # sanity check: the unreduced path with all spreads equal to COMMON must reproduce the common optimum
        sanity_cells = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.5)]
        jobs += [("sanity_common_unreduced", a, g, variant_spreads("common", unit)) for a, g in sanity_cells]
        jobs += [("sanity_common_reduced", a, g, None) for a, g in sanity_cells]
        phase1 = pool.map(solve, jobs, chunksize=1)
        by = {(r["tag"], r["alpha"], r["gamma"]): r for r in phase1}

        rows, differing = [], []
        for a, g in cells:
            c, u = by[("common", a, g)], by[("unit_specific", a, g)]
            diff = None if u["status"] != "completed" else (u["cost"] != c["cost"] or u["chromosome"] != c["chromosome"])
            row = {"alpha": a, "gamma": g, "common": c, "unit_specific": u, "diff": diff}
            if u["status"] == "completed":
                ts = library.get(tuple(u["chromosome"]), [])
                u["in_crisp_library"] = bool(ts)
                u["matching_T_range"] = [min(ts), max(ts)] if ts else None
            rows.append(row)
            if diff:
                differing.append((a, g))

        print(f"Phase 2: decomposition and sensitivity in {len(differing)} differing cells", flush=True)
        jobs = []
        for a, g in differing:
            for v in ("thermal_only", "hydro_only"):
                jobs.append((f"decomp:{v}", a, g, variant_spreads(v, unit)))
            for dh in HYDRO_SWEEP:
                s = dict(unit); s[HYDRO_MW] = (dh, 0.0)
                jobs.append((f"hydro:{dh}", a, g, s))
            for name, ct_l, ct_r, st_l, st_r in THERMAL_PERTURBATIONS:
                s = dict(unit)
                for mw in s:
                    if mw == 20:
                        s[mw] = (ct_l, ct_r)
                    elif mw != HYDRO_MW:
                        s[mw] = (st_l, st_r)
                jobs.append((f"thermal:{name}", a, g, s))
        phase2 = pool.map(solve, jobs, chunksize=1)

    p2 = {}
    for r in phase2:
        p2.setdefault((r["alpha"], r["gamma"]), []).append(r)
    for row in rows:
        key = (row["alpha"], row["gamma"])
        if key not in p2:
            continue
        rs = p2[key]
        row["decomposition"] = {
            "common": row["common"]["cost"],
            "thermal_only": next(r["cost"] for r in rs if r["tag"] == "decomp:thermal_only"),
            "hydro_only": next(r["cost"] for r in rs if r["tag"] == "decomp:hydro_only"),
            "unit_specific": row["unit_specific"]["cost"],
        }
        row["sensitivity"] = {
            "hydro": [{"delta_hydro": float(r["tag"].split(":")[1]), "status": r["status"],
                       "cost": r["cost"], "fleet": r["fleet"]} for r in rs if r["tag"].startswith("hydro:")],
            "thermal": [{"perturbation": r["tag"].split(":", 1)[1], "status": r["status"],
                         "cost": r["cost"], "fleet": r["fleet"]} for r in rs if r["tag"].startswith("thermal:")],
        }

    sanity = [{"alpha": a, "gamma": g,
               "reduced_cost": by[("sanity_common_reduced", a, g)]["cost"],
               "unreduced_cost": by[("sanity_common_unreduced", a, g)]["cost"],
               "agree": by[("sanity_common_reduced", a, g)]["chromosome"] == by[("sanity_common_unreduced", a, g)]["chromosome"]}
              for a, g in sanity_cells]

    completed = [r for r in rows if r["unit_specific"]["status"] == "completed"]
    times = [r["unit_specific"]["runtime_seconds"] for r in completed]
    ctimes = [r["common"]["runtime_seconds"] for r in rows]
    summary = {
        "n_cells": len(rows),
        "n_completed": len(completed),
        "n_timeout": len(rows) - len(completed),
        "n_differing": sum(1 for r in rows if r["diff"]),
        "n_differing_in_crisp_library": sum(1 for r in rows if r["diff"] and r["unit_specific"].get("in_crisp_library")),
        "unreduced_runtime_seconds": [min(times), max(times)] if times else None,
        "common_runtime_seconds": [min(ctimes), max(ctimes)],
        "crisp_library": {"n_thresholds": n_thresholds, "n_distinct_designs": len(library)},
        "time_limit_seconds": TIME_LIMIT,
        "sanity_unreduced_equals_reduced": sanity,
    }
    return {"description": __doc__.split("\n\n")[0],
            "assumptions": {"common_spread": COMMON, "unit_spreads": {str(k): v for k, v in unit.items()},
                            "source_of_unit_spreads": "hydro: RTS-79 Table 6 (90% available capacity in quarters 3-4); thermal: assumption of this study"},
            "summary": summary, "cells": rows}


def fmt(cost):
    return f"\\${cost}\\,M" if cost is not None else "timeout"


def table(report: dict) -> str:
    s = report["summary"]
    diff = [r for r in report["cells"] if r["diff"]]
    lines = [
        r"\begin{table}[!htbp]",
        r"\providecommand{\revtablecolor}{}\revtablecolor",
        r"\centering",
        r"\caption{Policy cells of the IEEE RTS-79 fleet in which unit-specific capacity spreads change the "
        r"optimal fleet relative to the common spread $\delta=0.10$. Unit-specific spreads: hydro "
        r"$(\delta_l,\delta_r)=(0.10,0)$ from RTS-79 Table~6; combustion turbine $(0.15,0.05)$ and steam units "
        r"$(0.05,0.02)$, an assumption of this study. ``Thermal only'' changes the thermal units and keeps hydro at "
        r"the common spread; ``hydro only'' "
        r"does the reverse.}",
        r"\label{tab:rts_unit_specific}",
        r"\footnotesize",
        r"\begin{tabular}{ccp{3.0cm}p{3.0cm}cccc}",
        r"\hline",
        r"$\alpha$ & $\gamma$ & Common-spread fleet & Unit-specific fleet & Common & Thermal only & Hydro only & Unit-specific\\",
        r"\hline",
    ]
    for r in diff:
        d = r["decomposition"]
        lines.append(f"{r['alpha']:.1f} & {r['gamma']:.1f} & {r['common']['fleet']} & {r['unit_specific']['fleet']} & "
                     f"{fmt(d['common'])} & {fmt(d['thermal_only'])} & {fmt(d['hydro_only'])} & {fmt(d['unit_specific'])}\\\\")
    lo, hi = s["unreduced_runtime_seconds"]
    note = (f"{s['n_differing']} of {s['n_cells']} cells differ; {s['n_timeout']} timed out at "
            f"{s['time_limit_seconds']:.0f}\\,s. In {s['n_differing_in_crisp_library']} of the {s['n_differing']} differing cells "
            f"the unit-specific optimum coincides with the crisp optimum of some integer threshold (library: "
            f"{s['crisp_library']['n_thresholds']} thresholds, {s['crisp_library']['n_distinct_designs']} distinct designs). "
            f"Unreduced solve time per cell: "
            f"${lo:.3f}$--${hi:.0f}$\\,s; reduced formulation: $\\le 0.01$\\,s per cell.")
    lines += [r"\hline",
              r"\multicolumn{8}{p{\dimexpr\linewidth-2\tabcolsep\relax}}{\textit{Note}: " + note + r"}\\",
              r"\hline", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    workers = int(os.environ.get("RTS_WORKERS", "2"))
    report = run(workers)
    (ROOT / "results" / "generated" / "rts_unit_specific_v7.json").write_text(json.dumps(report, indent=2) + "\n")
    tex = table(report)
    (ROOT / "results" / "generated" / "table_rts_unit_specific.tex").write_text(tex)
    (ROOT / "manuscript" / "table_rts_unit_specific.tex").write_text(tex)
    print(json.dumps(report["summary"], indent=2))
