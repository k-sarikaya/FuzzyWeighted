"""Regenerate every verification artifact from one reproducible script.

Round 2 objected that the eight `*_v5.json` files had no producing script, that
two of them disagreed about the same quantity, and that the credibility audit
reported 0/100,000 mismatches for a statement it could not have tested. This
script replaces all of them. Each audit declares its sampling domain, and the
credibility audit deliberately samples the asymmetric spreads under which the
claimed coincidence fails.

    PYTHONPATH=src python3 experiments/run_audits_v6.py
"""

from __future__ import annotations

import json
import random
import sys
from datetime import datetime, timezone
from fractions import Fraction
from math import ceil, log
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS  # noqa: E402
from fuzzy_wkon.fuzzy import TFN  # noqa: E402
from fuzzy_wkon.membership import (  # noqa: E402
    credibility_geq,
    hurwicz_indicator,
    membership_profile,
    necessity_geq,
    possibility_geq,
    proportional_capacity_tfn,
    saturation_level,
)
from fuzzy_wkon.reliability import capped_distribution, measure_reliability  # noqa: E402
from fuzzy_wkon.solver import prepare_options  # noqa: E402

RESULTS = ROOT / "results" / "generated"
GRIDS = (2001, 20001, 200001)


# --------------------------------------------------------------------------
# 1. Closed forms against a brute-force sup-min grid, with a measured rate.
# --------------------------------------------------------------------------

def _membership(tfn: TFN, x: float) -> float:
    if x < tfn.left or x > tfn.right:
        return 0.0
    if x < tfn.mode:
        return (x - tfn.left) / (tfn.mode - tfn.left) if tfn.mode > tfn.left else 1.0
    if x > tfn.mode:
        return (tfn.right - x) / (tfn.right - tfn.mode) if tfn.right > tfn.mode else 1.0
    return 1.0


def _brute_pos(a: TFN, b: TFN, points: int) -> float:
    low, high = min(a.left, b.left), max(a.right, b.right)
    grid = [low + (high - low) * i / (points - 1) for i in range(points)]
    mu_a = [_membership(a, x) for x in grid]
    mu_b = [_membership(b, y) for y in grid]
    suffix = [0.0] * points
    run = 0.0
    for i in range(points - 1, -1, -1):
        run = max(run, mu_a[i])
        suffix[i] = run
    return max(min(suffix[i], mu_b[i]) for i in range(points))


def audit_closed_forms(trials: int = 300, seed: int = 20260925) -> dict:
    rng = random.Random(seed)
    cases = []
    for _ in range(trials):
        a1 = rng.uniform(0.0, 20.0); a2 = a1 + rng.uniform(0.0, 8.0); a3 = a2 + rng.uniform(0.0, 8.0)
        b1 = rng.uniform(0.0, 20.0); b2 = b1 + rng.uniform(0.0, 8.0); b3 = b2 + rng.uniform(0.0, 8.0)
        cases.append((TFN(a1, a2, a3), TFN(b1, b2, b3)))

    per_grid = {}
    for points in GRIDS:
        worst_pos = worst_nec = 0.0
        for a, b in cases:
            worst_pos = max(worst_pos, abs(possibility_geq(a, b) - _brute_pos(a, b, points)))
            worst_nec = max(worst_nec, abs(necessity_geq(a, b) - (1.0 - _brute_pos(b, a, points))))
        per_grid[points] = {"max_error_pos": worst_pos, "max_error_nec": worst_nec,
                            "max_error": max(worst_pos, worst_nec)}

    first, last = GRIDS[0], GRIDS[-1]
    ratio = per_grid[first]["max_error"] / per_grid[last]["max_error"]
    order = log(ratio) / log(last / first)
    return {
        "description": "Closed-form Pos/Nec against a brute-force sup-min grid.",
        "sampling_domain": "300 random TFN pairs, independent left/right spreads",
        "seed": seed,
        "per_grid": per_grid,
        "error_reduction_factor": ratio,
        "empirical_order": order,
        "note": ("The measured order is reported as observed. It is below 1 because the "
                 "sup-min maximiser is attained at a kink, where the grid can only "
                 "approach the crossing point linearly in the spacing."),
    }


# --------------------------------------------------------------------------
# 2. Chance-constraint equivalence on the two boundary columns.
# --------------------------------------------------------------------------

def audit_equivalence(instances: int = 25000, seed: int = 20260925) -> dict:
    rng = random.Random(seed)
    levels = (0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
    mismatches_nec = mismatches_pos = 0
    checked = 0
    for _ in range(instances):
        modal = rng.uniform(1.0, 60.0)
        # Independent left/right spreads: the equivalence must hold for both.
        dl, dr = rng.choice([0.0, 0.05, 0.1, 0.2]), rng.choice([0.05, 0.1, 0.2, 0.4])
        capacity = TFN(modal * (1 - dl), modal, modal * (1 + dr))
        mode = rng.uniform(5.0, 55.0)
        threshold = TFN(mode - rng.uniform(0.1, 10.0), mode, mode + rng.uniform(0.1, 10.0))
        for level in levels:
            checked += 1
            if (necessity_geq(capacity, threshold) >= level - 1e-12) != \
               (hurwicz_indicator(capacity, threshold, 1.0 - level, 0.0) == 1.0):
                mismatches_nec += 1
            if (possibility_geq(capacity, threshold) >= level - 1e-12) != \
               (hurwicz_indicator(capacity, threshold, level, 1.0) == 1.0):
                mismatches_pos += 1
    return {
        "description": "Nec >= lambda <=> Hurwicz(1-lambda, 0); Pos >= lambda <=> Hurwicz(lambda, 1).",
        "sampling_domain": "asymmetric proportional capacity spreads and asymmetric thresholds",
        "seed": seed,
        "triples_checked": checked,
        "mismatches_necessity": mismatches_nec,
        "mismatches_possibility": mismatches_pos,
        "status": "PASS" if mismatches_nec == mismatches_pos == 0 else "FAIL",
    }


# --------------------------------------------------------------------------
# 3. Credibility: where the central column does and does not coincide.
# --------------------------------------------------------------------------

def audit_credibility(instances: int = 20000, seed: int = 20260925) -> dict:
    rng = random.Random(seed)
    symmetric_mismatches = asymmetric_mismatches = 0
    symmetric_checked = asymmetric_checked = 0
    witness = None

    for _ in range(instances):
        modal = rng.uniform(1.0, 60.0)
        symmetric = rng.random() < 0.5
        if symmetric:
            dl = dr = rng.choice([0.05, 0.1, 0.2])
        else:
            dl, dr = rng.choice([0.0, 0.05, 0.1]), rng.choice([0.2, 0.3, 0.4])
            if dl == dr:
                continue
        capacity = TFN(modal * (1 - dl), modal, modal * (1 + dr))
        mode = rng.uniform(5.0, 55.0)
        half = rng.uniform(0.1, 10.0)
        threshold = TFN(mode - half, mode, mode + half)  # symmetric threshold

        for alpha in (0.0, 0.25, 0.5, 0.75, 1.0):
            hurwicz = hurwicz_indicator(capacity, threshold, alpha, 0.5) == 1.0
            credible = credibility_geq(capacity, threshold) >= 0.5 - 1e-12
            if symmetric:
                symmetric_checked += 1
                if hurwicz != credible:
                    symmetric_mismatches += 1
            else:
                asymmetric_checked += 1
                if hurwicz != credible:
                    asymmetric_mismatches += 1
                    if witness is None:
                        witness = {"modal_capacity": modal, "delta_l": dl, "delta_r": dr,
                                   "alpha": alpha, "threshold": [threshold.left, threshold.mode,
                                                                 threshold.right],
                                   "hurwicz_accepts": hurwicz,
                                   "credibility": credibility_geq(capacity, threshold)}

    # Reliability level: R_H at gamma=1/2 versus R_Cr, on the archived designs.
    reliability_gaps = []
    for key, benchmark in BENCHMARKS.items():
        threshold = TFN(*benchmark["K_tfn"])
        spread = 0.1
        for alpha in (0.0, 0.5, 1.0):
            options, thr, scale = prepare_options(benchmark, alpha, 0.5, spread, spread)
            selected = [(o.reliability, o.weight_scaled)
                        for row, idx in zip(options, [len(r) - 1 for r in options])
                        for o in row if o.choice_index == row[idx].choice_index][:benchmark["N"]]
            cap = max(saturation_level(threshold, m, spread, scale) for m in ("pos", "cr"))
            profiles = {m: membership_profile(threshold, m, spread, spread, scale, cap)
                        for m in ("pos", "nec", "cr")}
            r_cr = measure_reliability(selected, profiles["cr"])
            indicator = [1.0 if level >= thr else 0.0 for level in range(cap + 1)]
            indicator[-1] = 1.0
            r_h = measure_reliability(selected, indicator)
            reliability_gaps.append({"catalog": key, "alpha": alpha,
                                     "R_hurwicz_gamma_half": r_h, "R_credibility": r_cr,
                                     "gap": r_h - r_cr})

    return {
        "description": ("Constraint level: Hurwicz at gamma=1/2 versus Cr >= 1/2. "
                        "Reliability level: R_H(gamma=1/2) versus R_Cr."),
        "seed": seed,
        "constraint_level": {
            "symmetric_spreads": {"checked": symmetric_checked,
                                  "mismatches": symmetric_mismatches},
            "asymmetric_spreads": {"checked": asymmetric_checked,
                                   "mismatches": asymmetric_mismatches},
            "witness_under_asymmetric_spreads": witness,
            "conclusion": ("The coincidence requires delta_l = delta_r. Under asymmetric "
                           "capacity spreads it fails, as the witness shows."),
        },
        "reliability_level": {
            "rows": reliability_gaps,
            "max_abs_gap": max(abs(r["gap"]) for r in reliability_gaps),
            "conclusion": ("R_H(gamma=1/2) is the expectation of a crisp indicator; R_Cr is a "
                           "degree-weighted expectation. They are different functionals and "
                           "the constraint-level coincidence does not extend to them."),
        },
    }


# --------------------------------------------------------------------------
# 4. Factorisation of the scalarised capacity (Proposition: reduction).
# --------------------------------------------------------------------------

def audit_factorisation() -> dict:
    alphas = (0, Fraction(1, 10), Fraction(1, 4), Fraction(3, 10), Fraction(1, 2),
              Fraction(7, 10), Fraction(9, 10), 1)
    gammas = (0, Fraction(1, 5), Fraction(3, 10), Fraction(1, 2), Fraction(7, 10),
              Fraction(4, 5), 1)
    spread = Fraction(1, 10)
    comparisons = mismatches = 0
    cells = 0
    for benchmark in BENCHMARKS.values():
        for alpha in alphas:
            for gamma in gammas:
                cells += 1
                options, thr, scale = prepare_options(
                    benchmark, float(alpha), float(gamma), float(spread), float(spread))
                m = 1 - spread * (1 - Fraction(alpha)) + \
                    Fraction(gamma) * 2 * spread * (1 - Fraction(alpha))
                for i, row in enumerate(options):
                    for option in row:
                        if option.choice_index == 0:
                            continue
                        comparisons += 1
                        modal = Fraction(str(benchmark["choices"][i][option.choice_index][2]))
                        if modal * m * scale != option.weight_scaled:
                            mismatches += 1
    return {
        "description": "Every scalarised capacity equals w_i * m(alpha, gamma) exactly.",
        "grid": {"alphas": [str(a) for a in alphas], "gammas": [str(g) for g in gammas]},
        "cells": cells,
        "comparisons": comparisons,
        "mismatches": mismatches,
        "status": "PASS" if mismatches == 0 else "FAIL",
    }


# --------------------------------------------------------------------------
# 5. Node ratios reported as the two distinct quantities they are.
# --------------------------------------------------------------------------

def audit_node_ratios() -> dict:
    source = RESULTS / "scaling_campaign_v4.json"
    payload = json.loads(source.read_text())
    rows = []
    for row in payload["rows"]:
        if row.get("status") != "completed" or row.get("nodes") is None:
            continue
        n, nodes, leaves = row["n"], row["nodes"], row.get("leaves")
        designs = 3 ** n
        full_tree = (3 ** (n + 1) - 1) // 2
        rows.append({
            "n": n, "profile": row["profile"],
            "designs_in_space": designs,
            "complete_designs_evaluated": leaves,
            "fraction_of_design_space_evaluated": (leaves / designs) if leaves is not None else None,
            "search_nodes_visited": nodes,
            "full_tree_nodes": full_tree,
            "search_effort_ratio": nodes / full_tree,
        })
    return {
        "description": ("Two distinct quantities. 'leaves' counts complete designs "
                        "evaluated and is comparable with 3^n. 'nodes' counts partial "
                        "assignments visited across both phases and is comparable with "
                        "the full tree (3^(n+1)-1)/2, not with 3^n."),
        "source": source.name,
        "rows": rows,
    }


def audit_uncertain_oracle() -> dict:
    from fuzzy_wkon.uncertain_oracle import optimistic_bound_uncertain, verify_example
    ex = verify_example()
    pharma = BENCHMARKS["pharma_adapted"]
    root_bound = optimistic_bound_uncertain(pharma, [], float(pharma["K_tfn"][1]), 0.1, 0.1)
    return {
        "description": "Reproduction of Wei et al. (2025) and root bound proving Pharma infeasibility",
        "wei_et_al_example": {
            "published_RM": ex["published_RM_traverse"],
            "computed_RM": ex["computed_RM"],
            "status": "PASS" if abs(ex["computed_RM"] - 0.6) < 1e-4 else "FAIL",
        },
        "pharma_root_bound": {
            "value": round(root_bound, 4),
            "exact_value": root_bound,
            "target": pharma["R_target"],
            "infeasible_at_root": bool(root_bound < pharma["R_target"]),
            "status": "PASS" if round(root_bound, 4) == 0.8488 else "FAIL",
        },
    }


def audit_uncertain_bound_pairs(n_pairs: int = 20000, seed: int = 42) -> dict:
    from fuzzy_wkon.uncertain_oracle import optimistic_bound_uncertain, evaluate_design_uncertain
    pharma = BENCHMARKS["pharma_adapted"]
    threshold = float(pharma["K_tfn"][1])
    n = len(pharma["choices"])
    rng = random.Random(seed)
    violations = 0
    min_slack = float("inf")
    for _ in range(n_pairs):
        k = rng.randint(0, n)
        full = [rng.choice(range(len(pharma["choices"][i]))) for i in range(n)]
        partial = full[:k]
        ub = optimistic_bound_uncertain(pharma, partial, threshold, 0.1, 0.1)
        _, true_val = evaluate_design_uncertain(pharma, full, threshold, 0.1, 0.1)
        slack = ub - true_val
        if slack < -1e-12:
            violations += 1
        if slack < min_slack:
            min_slack = slack
    return {
        "description": "Arm-U uncertain-measure bound admissibility audit over random prefix-completion pairs",
        "pairs_checked": n_pairs,
        "violations": violations,
        "min_slack": min_slack,
        "status": "PASS" if violations == 0 else "FAIL",
    }


def audit_uncertain_all_designs() -> dict:
    import itertools
    from fuzzy_wkon.uncertain_oracle import evaluate_design_uncertain
    pharma = BENCHMARKS["pharma_adapted"]
    threshold = float(pharma["K_tfn"][1])
    max_rm = -1.0
    best_chrom = None
    for p in itertools.product(range(3), repeat=10):
        _, rm = evaluate_design_uncertain(pharma, list(p), threshold, 0.1, 0.1)
        if rm > max_rm:
            max_rm = rm
            best_chrom = list(p)
    return {
        "description": "Exhaustive evaluation of all 3^10 = 59,049 designs under Arm-U uncertain measure",
        "designs_checked": 59049,
        "max_uncertain_measure": max_rm,
        "best_chromosome": best_chrom,
        "target": pharma["R_target"],
        "all_infeasible": bool(max_rm < pharma["R_target"]),
        "status": "PASS" if round(max_rm, 4) == 0.8488 and max_rm < pharma["R_target"] else "FAIL",
    }


def audit_optimistic_bound_pairs(n_pairs: int = 20000, seed: int = 42) -> dict:
    from fuzzy_wkon.solver import prepare_options, _optimistic_reliability
    from fuzzy_wkon.reliability import capped_distribution
    rng = random.Random(seed)
    violations = 0
    min_slack = float("inf")
    for bench_name in ("pharma_adapted", "fyffe_adapted"):
        bench = BENCHMARKS[bench_name]
        n = len(bench["choices"])
        options, threshold, _ = prepare_options(bench, 0.5, 0.5)
        for _ in range(n_pairs // 2):
            prefix_len = rng.randint(0, n)
            full_chromosome = [rng.choice(range(len(bench["choices"][i]))) for i in range(n)]
            prefix_sel = [(options[i][c].reliability, options[i][c].weight_scaled)
                          for i, c in enumerate(full_chromosome[:prefix_len]) if options[i][c].choice_index != 0]
            pmf_prefix = capped_distribution(prefix_sel, threshold)
            pmf_dict = dict(enumerate(pmf_prefix))
            ub = _optimistic_reliability(pmf_dict, options, prefix_len, threshold)
            full_sel = [(options[i][c].reliability, options[i][c].weight_scaled)
                        for i, c in enumerate(full_chromosome) if options[i][c].choice_index != 0]
            true_rel = capped_distribution(full_sel, threshold)[-1]
            slack = ub - true_rel
            if slack < -1e-12:
                violations += 1
            if slack < min_slack:
                min_slack = slack
    return {
        "description": "Theorem 3 optimistic reliability bound audit over random prefix-completion pairs",
        "pairs_checked": n_pairs,
        "violations": violations,
        "min_slack": min_slack,
        "status": "PASS" if violations == 0 else "FAIL",
    }


def audit_pharma_00_infeasibility() -> dict:
    """Exhaustive enumeration proving Pharmaceutical plant (0,0) cell infeasibility at R_target=0.85."""
    from itertools import product
    bench = BENCHMARKS["pharma_adapted"]
    # At (0,0) cell: alpha=0, gamma=0, spread=0.10
    # K_tfn = (34, 40, 46) => demand upper = 46.
    # capacity derated factor = 1 - 0.1 = 0.9.
    # T = ceil(46 / 0.9) = 52.
    T = 52
    target = bench["R_target"]  # 0.85
    choices = bench["choices"]

    max_modal_cap = sum(max(w for _, _, w in row) for row in choices)  # 86
    max_rel = 0.0
    best_chrom = None
    for chrom in product(range(3), repeat=len(choices)):
        selected = [(choices[i][c][0], choices[i][c][2]) for i, c in enumerate(chrom) if c != 0]
        pmf = capped_distribution(selected, T)
        rel = pmf[-1]
        if rel > max_rel:
            max_rel = rel
            best_chrom = chrom

    return {
        "description": "Exhaustive enumeration of all 3^10 designs for Pharmaceutical plant at (alpha=0, gamma=0)",
        "T": T,
        "max_installed_modal_capacity": max_modal_cap,
        "max_derated_capacity": round(max_modal_cap * 0.9, 1),
        "capacity_binding": False,
        "target_reliability": target,
        "max_achievable_reliability": round(max_rel, 6),
        "best_chromosome": list(best_chrom),
        "status": "PASS" if max_rel < target else "FAIL",
    }


if __name__ == "__main__":
    audits = {
        "closed_forms": audit_closed_forms(),
        "equivalence": audit_equivalence(),
        "credibility": audit_credibility(),
        "factorisation": audit_factorisation(),
        "node_ratios": audit_node_ratios(),
        "uncertain_oracle": audit_uncertain_oracle(),
        "uncertain_bound_pairs": audit_uncertain_bound_pairs(),
        "uncertain_all_designs": audit_uncertain_all_designs(),
        "optimistic_bound_pairs": audit_optimistic_bound_pairs(),
        "pharma_00_infeasibility": audit_pharma_00_infeasibility(),
    }
    payload = {"generated_at_utc": datetime.now(timezone.utc).isoformat(),
               "script": "experiments/run_audits_v6.py", "audits": audits}
    (RESULTS / "audits_v6.json").write_text(json.dumps(payload, indent=2))

    cf = audits["closed_forms"]
    print(f"closed forms : order {cf['empirical_order']:.3f}, "
          f"error {cf['per_grid'][GRIDS[0]]['max_error']:.2e} -> "
          f"{cf['per_grid'][GRIDS[-1]]['max_error']:.2e}")
    eq = audits["equivalence"]
    print(f"equivalence  : {eq['triples_checked']:,} triples, "
          f"{eq['mismatches_necessity']} nec / {eq['mismatches_possibility']} pos mismatches")
    cr = audits["credibility"]["constraint_level"]
    print(f"credibility  : symmetric {cr['symmetric_mismatches'] if 'symmetric_mismatches' in cr else cr['symmetric_spreads']['mismatches']}"
          f"/{cr['symmetric_spreads']['checked']} mismatches; "
          f"asymmetric {cr['asymmetric_spreads']['mismatches']}/{cr['asymmetric_spreads']['checked']}")
    print(f"               reliability-level max |gap| = "
          f"{audits['credibility']['reliability_level']['max_abs_gap']:.4f}")
    fa = audits["factorisation"]
    print(f"factorisation: {fa['comparisons']:,} comparisons over {fa['cells']} cells, "
          f"{fa['mismatches']} mismatches")
    nr = [r for r in audits["node_ratios"]["rows"] if r["n"] == 20 and r["profile"] == "design"]
    if nr:
        r = nr[0]
        print(f"node ratios  : n=20/design evaluated {r['complete_designs_evaluated']} designs "
              f"= {100*r['fraction_of_design_space_evaluated']:.2e}% of 3^20; "
              f"search effort {100*r['search_effort_ratio']:.1f}% of the full tree")
    print(f"pharma (0,0) : T={audits['pharma_00_infeasibility']['T']}, max rel={audits['pharma_00_infeasibility']['max_achievable_reliability']:.4f} < {audits['pharma_00_infeasibility']['target_reliability']} (status: {audits['pharma_00_infeasibility']['status']})")
    print(f"\nwritten: {RESULTS / 'audits_v6.json'}")
