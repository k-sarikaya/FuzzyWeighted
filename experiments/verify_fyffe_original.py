"""Reproduce the numerical Fyffe et al. (1968) Table III allocation."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.fyffe_original import SOURCE_TABLE, solve_original_fyffe


def main() -> None:
    result = solve_original_fyffe(cost_limit=130, weight_limit=170)
    published = {"cost": 119, "weight": 170, "reliability_reported": 0.9700}
    result["source"] = {
        "citation": "Fyffe, Hines and Lee (1968), Table I p. 68 and Table III p. 68",
        "table_csv": str(SOURCE_TABLE.relative_to(ROOT)),
        "table_csv_sha256": hashlib.sha256(SOURCE_TABLE.read_bytes().replace(b"\r\n", b"\n")).hexdigest()  # line-ending independent,
    }
    result["published_table_iii_check"] = {
        **published,
        "cost_matches": result["best_cost"] == published["cost"],
        "weight_matches": result["best_weight"] == published["weight"],
        "reliability_matches_to_4_decimals": round(result["best_reliability"], 4)
        == published["reliability_reported"],
    }
    result["claim_boundary"] = (
        "This verifies the original crisp series reliability-allocation example only. "
        "It neither validates nor changes the V2 adapted dual-fuzzy weighted catalog."
    )
    destination = ROOT / "results" / "generated" / "fyffe_original_reproduction.json"
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not all(result["published_table_iii_check"].values()):
        raise AssertionError("Original Fyffe reproduction does not match the published Table III")
    print(
        "Reproduced Fyffe Table III: "
        f"R={result['best_reliability']:.10f}, C={result['best_cost']}, W={result['best_weight']}"
    )


if __name__ == "__main__":
    main()
