"""Optimization trace of Supplementary Material S1 (three-subsystem example).

Follows Algorithm 1 of the main article as stated there: Phase 0 (effective values and suffix
tables), the greedy incumbent seed, Phase 1 with both pruning tests applied at the point of branching,
and Phase 2 on the iso-cost set.  Capacities and the threshold are exact rationals; reliabilities are
exact rationals as well, so no tolerance is needed.  This is an independent transcription of the
algorithm for the worked example, not the production solver (which requires a null option in every
catalog).

Writes results/generated/supp_trace_v11.json.

    python3 experiments/supp_trace_v11.py
"""
from __future__ import annotations

import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "generated" / "supp_trace_v11.json"

ALPHA, GAMMA = F(1, 2), F(0)
K_TFN = (F(6), F(8), F(10))
R_TARGET = F(9, 10)
# catalogs: subsystem -> list of (label, cost, reliability, capacity TFN)
CATALOG = [
    [(1, 10, F(95, 100), (F("4.32"), F("4.80"), F("5.28"))),
     (2, 20, F(97, 100), (F("4.95"), F("5.50"), F("6.05")))],
    [(1, 15, F(90, 100), (F("4.50"), F("5.00"), F("5.50")))],
    [(1, 15, F(85, 100), (F("4.50"), F("5.00"), F("5.50")))],
]


def cut(tfn, a):
    a1, a2, a3 = tfn
    return a1 + a * (a2 - a1), a3 - a * (a3 - a2)


def hurwicz_w(tfn):
    lo, up = cut(tfn, ALPHA)
    return GAMMA * up + (1 - GAMMA) * lo


def hurwicz_k():
    lo, up = cut(K_TFN, ALPHA)
    return up - GAMMA * (up - lo)


K = hurwicz_k()
OPTS = [[{"label": lab, "cost": c, "p": p, "w": hurwicz_w(t)} for lab, c, p, t in row] for row in CATALOG]
N = len(OPTS)
W_MAX_SUFFIX = [sum(max(o["w"] for o in OPTS[l]) for l in range(i + 1, N)) for i in range(N)]


def z_suf(i, q):
    """Minimum cost of subsystems i+1..n among completions supplying at least q (0-based i)."""
    rest = OPTS[i + 1:]
    best = None
    for combo in product(*rest):
        if sum(o["w"] for o in combo) >= q:
            c = sum(o["cost"] for o in combo)
            best = c if best is None else min(best, c)
    return best  # None = no completion supplies q


def convolve(pmf, p, w):
    nxt = {}
    for cap, mass in pmf.items():
        nxt[cap] = nxt.get(cap, 0) + mass * (1 - p)
        hi = min(K, cap + w)
        nxt[hi] = nxt.get(hi, 0) + mass * p
    return nxt


def rel(pmf):
    return pmf.get(K, F(0))


def optimistic(pmf, i):
    for l in range(i + 1, N):
        pmf = convolve(pmf, max(o["p"] for o in OPTS[l]), max(o["w"] for o in OPTS[l]))
    return rel(pmf)


def design_rel(design):
    pmf = {F(0): F(1)}
    for o in design:
        pmf = convolve(pmf, o["p"], o["w"])
    return rel(pmf)


def greedy_seed(log):
    design = [max(row, key=lambda o: o["w"]) for row in OPTS]
    r = design_rel(design)
    log.append({"step": "greedy start", "design": [o["label"] for o in design],
                "cost": sum(o["cost"] for o in design), "R": float(r)})
    if r < R_TARGET:
        return None
    improved = True
    while improved:
        improved = False
        for i, row in enumerate(OPTS):
            for o in sorted(row, key=lambda o: o["cost"]):
                if o["cost"] < design[i]["cost"]:
                    trial = design[:i] + [o] + design[i + 1:]
                    rt = design_rel(trial)
                    if rt >= R_TARGET:
                        design, improved = trial, True
                        log.append({"step": "greedy substitution", "design": [x["label"] for x in design],
                                    "cost": sum(x["cost"] for x in design), "R": float(rt)})
                        break
    return sum(o["cost"] for o in design)


def search(phase, log, best):
    def rec(i, z_part, w_part, pmf, x):
        if i == N:
            r = rel(pmf)
            log.append({"phase": phase, "event": "leaf", "design": x[:], "cost": z_part, "R": float(r)})
            if phase == 1 and r >= R_TARGET:
                best["Z"] = min(best["Z"], z_part)
            if phase == 2 and z_part == best["Z"] and r >= R_TARGET and r > best["R"]:
                best["R"], best["x"] = r, x[:]
            return
        for o in sorted(OPTS[i], key=lambda o: o["cost"]):
            z_next = z_part + o["cost"]
            q_next = max(F(0), K - w_part - o["w"])
            zs = z_suf(i, q_next)
            node = {"phase": phase, "node": x + [o["label"]], "Z_next": z_next, "q_next": float(q_next),
                    "Z_suf": zs, "W_bound": float(w_part + o["w"] + W_MAX_SUFFIX[i])}
            zc = best["Z"]
            if (zs is None or z_next + zs > zc) or (w_part + o["w"] + W_MAX_SUFFIX[i] < K):
                node["event"] = "pruned by Theorem 2"
                log.append(node)
                continue
            nxt = convolve(pmf, o["p"], o["w"])
            ub = optimistic(nxt, i)
            node["R_ub"] = float(ub)
            if ub < R_TARGET:
                node["event"] = "pruned by Theorem 3"
                log.append(node)
                continue
            node["event"] = "branch"
            log.append(node)
            x.append(o["label"])
            rec(i + 1, z_next, w_part + o["w"], nxt, x)
            x.pop()
    rec(0, 0, F(0), {F(0): F(1)}, [])


def main() -> dict:
    log: list[dict] = []
    seed = greedy_seed(log)
    best = {"Z": seed if seed is not None else float("inf"), "R": F(-1), "x": None}
    search(1, log, best)
    search(2, log, best)
    report = {
        "cell": [float(ALPHA), float(GAMMA)], "k_alpha_gamma": float(K), "R_target": float(R_TARGET),
        "effective_capacities": [[float(o["w"]) for o in row] for row in OPTS],
        "W_max_suffix": [float(v) for v in W_MAX_SUFFIX],
        "Z_suf_1": {str(float(q)): z_suf(0, q) for q in (F(0), K - OPTS[0][0]["w"], K - OPTS[0][1]["w"])},
        "R_designs": {"[1,1,1]": float(design_rel([OPTS[0][0], OPTS[1][0], OPTS[2][0]])),
                      "[2,1,1]": float(design_rel([OPTS[0][1], OPTS[1][0], OPTS[2][0]]))},
        "seed_cost": seed, "Z_star": best["Z"], "x_star": best["x"], "R_star": float(best["R"]),
        "log": log,
    }
    OUT.write_text(json.dumps(report, indent=1, default=str) + "\n")
    return report


if __name__ == "__main__":
    rep = main()
    for e in rep["log"]:
        print(e)
    print({k: v for k, v in rep.items() if k != "log"})
