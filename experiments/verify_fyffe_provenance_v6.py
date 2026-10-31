"""Provenance check for the Fyffe benchmark: transcription, reproduction, adaptation diff.

Three claims in the manuscript rest on this script.

1. The machine-readable archive of Table I of Fyffe, Hines and Lee (1968) matches
   the published table entry by entry. The reference values below were read from
   the published scan and are kept here so the check is reproducible without it.
2. Solving the source problem in its own model -- series architecture, integer
   parallel redundancy, cost and weight budgets -- recovers the published
   optimum of Table III.
3. The adapted catalog differs from the source in exactly the ways the manuscript
   states, and the count of verbatim rows is what it reports.

    PYTHONPATH=src python3 experiments/verify_fyffe_provenance_v6.py
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS                      # noqa: E402
from fuzzy_wkon.fyffe_original import load_table_i, solve_original_fyffe  # noqa: E402

# Table I as printed in IEEE Trans. Reliab. R-17(2), June 1968, p. 68.
# '*' in the source (no alternative available) appears here as a shorter row.
PUBLISHED_TABLE_I = {
    1:  [(0.90, 1, 3), (0.93, 1, 4), (0.91, 2, 2), (0.95, 2, 5)],
    2:  [(0.95, 2, 8), (0.94, 1, 10), (0.93, 1, 9)],
    3:  [(0.85, 2, 7), (0.90, 3, 5), (0.87, 1, 6), (0.92, 4, 4)],
    4:  [(0.83, 3, 5), (0.87, 4, 6), (0.85, 5, 4)],
    5:  [(0.94, 2, 4), (0.93, 2, 3), (0.95, 3, 5)],
    6:  [(0.99, 3, 5), (0.98, 3, 4), (0.97, 2, 5), (0.96, 2, 4)],
    7:  [(0.91, 4, 7), (0.92, 4, 8), (0.94, 5, 9)],
    8:  [(0.81, 3, 4), (0.90, 5, 7), (0.91, 6, 6)],
    9:  [(0.97, 2, 8), (0.99, 3, 9), (0.96, 4, 7), (0.91, 3, 8)],
    10: [(0.83, 4, 6), (0.85, 4, 5), (0.90, 5, 6)],
    11: [(0.94, 3, 5), (0.95, 4, 6), (0.96, 5, 6)],
    12: [(0.79, 2, 4), (0.82, 3, 5), (0.85, 4, 6), (0.90, 5, 7)],
    13: [(0.98, 2, 5), (0.99, 3, 5), (0.97, 2, 6)],
    14: [(0.90, 4, 6), (0.92, 4, 7), (0.95, 5, 6), (0.99, 6, 9)],
}

# Table III of the same paper: (design alternative, parallel units) per unit.
PUBLISHED_TABLE_III = {
    1: (3, 3), 2: (1, 2), 3: (4, 3), 4: (3, 3), 5: (2, 3), 6: (2, 2), 7: (1, 2),
    8: (1, 4), 9: (3, 2), 10: (2, 3), 11: (1, 2), 12: (1, 4), 13: (2, 2), 14: (3, 2),
}
PUBLISHED_OPTIMUM = {"cost": 119, "weight": 170, "reliability": 0.9700}


def check_transcription(table) -> dict:
    mismatches = []
    for unit, published in sorted(PUBLISHED_TABLE_I.items()):
        archived = [(a.reliability, a.cost, a.weight) for a in table[unit - 1]]
        if archived != published:
            mismatches.append({"unit": unit, "archived": archived, "published": published})
    return {"units": len(PUBLISHED_TABLE_I),
            "matching": len(PUBLISHED_TABLE_I) - len(mismatches),
            "mismatches": mismatches}


def check_reproduction() -> dict:
    out = solve_original_fyffe()
    policy = {d["functional_unit"]: (d["alternative"], d["parallel_units"])
              for d in out["design"]}
    disagreeing = {u: {"reproduced": policy[u], "published": PUBLISHED_TABLE_III[u]}
                   for u in PUBLISHED_TABLE_III if policy[u] != PUBLISHED_TABLE_III[u]}
    return {
        "published": PUBLISHED_OPTIMUM,
        "reproduced": {"cost": out["best_cost"], "weight": out["best_weight"],
                       "reliability": out["best_reliability"]},
        "cost_matches": out["best_cost"] == PUBLISHED_OPTIMUM["cost"],
        "weight_matches": out["best_weight"] == PUBLISHED_OPTIMUM["weight"],
        "reliability_matches_to_published_precision":
            round(out["best_reliability"], 4) == PUBLISHED_OPTIMUM["reliability"],
        "policy_units_disagreeing": disagreeing,
    }


def check_adaptation(table) -> dict:
    published = {u + 1: [(a.reliability, a.cost, a.weight) for a in alts]
                 for u, alts in enumerate(table)}
    verbatim, synthetic, null_options = 0, [], 0
    for unit, row in enumerate(BENCHMARKS["fyffe_adapted"]["choices"], start=1):
        for option in row:
            if tuple(option) == (0.0, 0, 0):
                null_options += 1
            elif tuple(option) in [tuple(x) for x in published[unit]]:
                verbatim += 1
            else:
                synthetic.append({"unit": unit, "option": list(option)})
    return {
        "non_null_options": verbatim + len(synthetic),
        "verbatim_from_table_i": verbatim,
        "synthetic_options": len(synthetic),
        "null_options_added": null_options,
        "synthetic_detail": synthetic,
        "structural_changes": [
            "the third field is reinterpreted from physical weight to operational capacity",
            "integer parallel redundancy is removed; one option per subsystem",
            "a (0,0,0) null option is added, which the source does not admit",
            "the source budgets C<=130 and W<=170 are not enforced",
        ],
    }


def main() -> None:
    table = load_table_i()
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "Fyffe, D. E., Hines, W. W., & Lee, N. K. (1968). "
                  "IEEE Transactions on Reliability, R-17(2), 64-69.",
        "transcription": check_transcription(table),
        "reproduction": check_reproduction(),
        "adaptation": check_adaptation(table),
    }
    dest = ROOT / "results" / "generated" / "fyffe_provenance_v6.json"
    dest.write_text(json.dumps(report, indent=2))

    t, r, a = report["transcription"], report["reproduction"], report["adaptation"]
    print(f"transcription: {t['matching']}/{t['units']} units match the published Table I")
    print(f"reproduction : C={r['reproduced']['cost']} W={r['reproduced']['weight']} "
          f"R={r['reproduced']['reliability']:.8f}  "
          f"(published C=119 W=170 R=0.9700); "
          f"{14 - len(r['policy_units_disagreeing'])}/14 units match Table III")
    print(f"adaptation   : {a['verbatim_from_table_i']} of {a['non_null_options']} "
          f"non-null options verbatim, {a['synthetic_options']} synthetic, "
          f"{a['null_options_added']} null options added")
    print("written:", dest)


if __name__ == "__main__":
    main()
