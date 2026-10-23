"""Generate V2 solver artifacts from the frozen model contract."""

from __future__ import annotations

import hashlib
import json
import sys
import argparse
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS
from fuzzy_wkon.solver import solve_exact


def config_signature(config: dict) -> str:
    payload = json.dumps(config, sort_keys=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", choices=sorted(BENCHMARKS), default=None,
                        help="run one adapted benchmark; omit for the full grid")
    parser.add_argument("--spread", default=None,
                        help="spread name from v2_default.json; default_spread otherwise")
    args = parser.parse_args()
    config_path = ROOT / "configs" / "v2_default.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    spread_by_name = {row["name"]: row for row in config["spreads"]}
    default = spread_by_name[args.spread or config["default_spread"]]
    output_dir = ROOT / "results" / "generated"
    output_dir.mkdir(parents=True, exist_ok=True)
    campaign = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_signature": config_signature(config),
        "objective": config["solver_objective"],
        "model": config["reliability_target_semantics"],
        "results": [],
    }
    selected_benchmarks = BENCHMARKS
    if args.benchmark is not None:
        selected_benchmarks = {args.benchmark: BENCHMARKS[args.benchmark]}
    for benchmark_name, benchmark in selected_benchmarks.items():
        for alpha in config["alpha_grid"]:
            for gamma in config["gamma_grid"]:
                result = solve_exact(benchmark, alpha, gamma, default["left"], default["right"])
                result.update({"benchmark": benchmark_name, "spread": default["name"]})
                campaign["results"].append(result)
                print(json.dumps({"benchmark": benchmark_name, "alpha": alpha, "gamma": gamma,
                                  "cost": result["best_cost"], "feasible": result["feasible"]}))
    output = output_dir / "campaign_v2_default.json"
    output.write_text(json.dumps(campaign, indent=2), encoding="utf-8")
    print(f"Wrote {len(campaign['results'])} cells to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
