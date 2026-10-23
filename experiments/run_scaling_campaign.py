"""Run the declared bounded scaling campaign and export a machine-readable audit."""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.scaling import PROFILES, synthetic_benchmark
from fuzzy_wkon.solver import solve_exact


def _solve(payload: tuple[int, str, int, float, float]) -> dict:
    n, profile, seed, alpha, gamma = payload
    benchmark = synthetic_benchmark(n, profile, seed)
    result = solve_exact(benchmark, alpha, gamma)
    return {
        "n": n,
        "profile": profile,
        "seed": seed,
        "alpha": alpha,
        "gamma": gamma,
        "status": "completed",
        "feasible": result["feasible"],
        "best_cost": result["best_cost"],
        "best_reliability": result["best_reliability"],
        "runtime_seconds": result["runtime_seconds"],
        "nodes": result["stats"]["nodes"],
        "leaves": result["stats"]["leaves"],
        "search_space": 3 ** n,
        "generator": benchmark["generator"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=float, default=12.0)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--sizes", nargs="*", type=int, default=[8, 10, 12, 16, 20, 24, 28, 32])
    args = parser.parse_args()
    alpha, gamma = 0.3, 0.7
    rows = []
    context = mp.get_context("spawn")
    for n in args.sizes:
        for profile in PROFILES:
            payload = (n, profile, args.seed, alpha, gamma)
            with context.Pool(1) as pool:
                pending = pool.apply_async(_solve, (payload,))
                try:
                    rows.append(pending.get(timeout=args.timeout))
                except mp.TimeoutError:
                    pool.terminate()
                    rows.append({
                        "n": n,
                        "profile": profile,
                        "seed": args.seed,
                        "alpha": alpha,
                        "gamma": gamma,
                        "status": "timeout",
                        "timeout_seconds": args.timeout,
                        "search_space": 3 ** n,
                        "generator": synthetic_benchmark(n, profile, args.seed)["generator"],
                    })
                except Exception as error:
                    pool.terminate()
                    rows.append({
                        "n": n,
                        "profile": profile,
                        "seed": args.seed,
                        "alpha": alpha,
                        "gamma": gamma,
                        "status": "error",
                        "error_type": type(error).__name__,
                        "error_message": str(error),
                        "search_space": 3 ** n,
                        "generator": synthetic_benchmark(n, profile, args.seed)["generator"],
                    })
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "bounded computational profile, not a general scalability claim",
        "solver": "two-phase exact branch-and-bound",
        "sizes": args.sizes,
        "profiles": PROFILES,
        "timeout_seconds_per_instance": args.timeout,
        "memory_measurement": "omitted because allocation tracing materially perturbed exact-search runtime",
        "runtime_environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "processor_environment": platform.processor() or "not exposed",
        },
        "rows": rows,
    }
    destination = ROOT / "results" / "generated" / "scaling_campaign.json"
    destination.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} scaling observations to {destination}")


if __name__ == "__main__":
    mp.freeze_support()
    main()
