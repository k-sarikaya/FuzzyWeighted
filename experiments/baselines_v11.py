"""Model baselines (v11): what a planner obtains by ignoring a source of imprecision, or by
treating it probabilistically, when the declared rule is the dual-fuzzy Hurwicz rule.

Decision-level baselines (all three instances; no new optimization)
--------------------------------------------------------------------
At a posture (alpha, gamma) the dual-fuzzy problem is the crisp problem at the integer threshold
T_d = ceil(kappa_d), kappa_d = k_H / m with m = 1 + d (1 - alpha)(2 gamma - 1) (Theorem 1).  A planner
who uses the same data and the same posture but a simpler model solves the crisp problem at

    crisp       T_b = ceil(k_2)                       (modal capacities and modal demand)
    w only      T_b = ceil(k_2 / m)                   (fuzzy capacities, crisp modal demand)
    k only      T_b = ceil(k_H)                       (crisp modal capacities, fuzzy demand)

Every T_b lies in the range of the parametric solution, so the baseline design is read from the
archived library.  It is then evaluated at T_d: it either fails the reliability target, or meets it
at a cost at least the optimum Z*(T_d).  At fixed alpha < 1 each kappa is monotone in gamma, so the
region {T_d = t, T_b = s} is a union of gamma intervals with closed-form ends; areas are integrated in
alpha by the midpoint rule.  Points at which no design meets the target at T_d are excluded and
reported separately.

Probabilistic counterpart (IEEE RTS-79 and Pharmaceutical-derived; exhaustive enumeration)
------------------------------------------------------------------------------------------
Unit states X_i are Bernoulli(p_i) as before; a common derating factor M and the demand K are random
and independent of X and of each other.  Success is M sum_i w_i X_i >= K, i.e. sum_i w_i X_i >= K / M,
and since the modal capacities are integers,

    R_prob(x) = sum_T q_T R(x, T),   q_T = P(ceil(K / M) = T),

where R(x, T) is the crisp reliability at threshold T.  The design minimizes cost subject to
R_prob(x) >= R_target, ties broken by the larger R_prob (the lexicographic rule of the main model).
Two distribution choices on the same supports: triangular (the membership functions read as densities)
and uniform.  The catalogs have 3^10 = 59,049 designs and are enumerated; the Fyffe-derived catalog
(4^14 designs) is not, because the exact solver's bounds are stated for a single threshold.

Writes results/generated/baselines_v11.json, manuscript/table_baselines.tex (main article) and
manuscript/table_baselines_supp.tex (Supplementary Material S7).

    python experiments/baselines_v11.py
"""
from __future__ import annotations

import itertools
import json
import sys
from math import ceil, floor
from pathlib import Path

import numpy as np
from scipy import integrate

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "generated"
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS          # noqa: E402
from fuzzy_wkon.reliability import static_reliability  # noqa: E402

SPREAD = 0.10
N_ALPHA = 20_000
EPS = 1e-12
INSTANCES = {"rts": "power_rts_adapted", "fyffe": "fyffe_adapted", "pharma": "pharma_adapted"}
PROBABILISTIC = ("rts", "pharma")


# ----------------------------------------------------------------------------- library and checks
def library(key: str) -> dict[int, dict]:
    if key == "rts":
        path = json.loads((RES / "parametric_stability_rts.json").read_text())["parametric_path"]
    else:
        path = json.loads((RES / "parametric_stability.json").read_text())["catalogs"][INSTANCES[key]][
            "parametric_path"]
    return {r["T"]: r for r in path}


def units(bench: dict, design) -> list[tuple[float, int]]:
    return [(bench["choices"][i][j][0], int(bench["choices"][i][j][2])) for i, j in enumerate(design) if j != 0]


def cost(bench: dict, design) -> int:
    return sum(bench["choices"][i][j][1] for i, j in enumerate(design))


def meets(bench: dict, design, T: int) -> bool:
    return static_reliability(units(bench, design), T) + EPS >= bench["R_target"]


# ----------------------------------------------------------------------------- effective thresholds
def k_hurwicz(k_tfn, alpha, gamma):
    k1, k2, k3 = k_tfn
    lo, up = k1 + alpha * (k2 - k1), k3 - alpha * (k3 - k2)
    return up - gamma * (up - lo)


def factor(alpha, gamma):
    return 1.0 + SPREAD * (1.0 - alpha) * (2.0 * gamma - 1.0)


def kappa(model: str, k_tfn, alpha, gamma):
    if model == "dual":
        return k_hurwicz(k_tfn, alpha, gamma) / factor(alpha, gamma)
    if model == "crisp":
        return float(k_tfn[1])
    if model == "w_only":
        return k_tfn[1] / factor(alpha, gamma)
    if model == "k_only":
        return k_hurwicz(k_tfn, alpha, gamma)
    raise ValueError(model)


def gamma_where(model: str, k_tfn, alpha: float, level: float) -> float | None:
    """The gamma at which kappa_model(alpha, gamma) = level (each kappa is monotone in gamma)."""
    k1, k2, k3 = k_tfn
    lo, up = k1 + alpha * (k2 - k1), k3 - alpha * (k3 - k2)
    d = SPREAD * (1.0 - alpha)
    if model == "dual":
        return (up - level * (1.0 - d)) / ((up - lo) + 2.0 * level * d)
    if model == "w_only":
        return ((k2 / level - 1.0) / d + 1.0) / 2.0 if d > 0 else None
    if model == "k_only":
        return (up - level) / (up - lo) if up > lo else None
    return None


def joint_areas(k_tfn, model: str, levels: range) -> dict[tuple[int, int], float]:
    """Area of {(alpha, gamma) : ceil(kappa_dual) = t, ceil(kappa_model) = s} for all (t, s)."""
    out: dict[tuple[int, int], float] = {}
    h = 1.0 / N_ALPHA
    for i in range(N_ALPHA):
        a = (i + 0.5) * h
        cuts = {0.0, 1.0}
        for L in levels:
            for m in ("dual", model):
                g = gamma_where(m, k_tfn, a, float(L))
                if g is not None and 0.0 < g < 1.0:
                    cuts.add(g)
        cuts = sorted(cuts)
        for g0, g1 in zip(cuts[:-1], cuts[1:]):
            gm = 0.5 * (g0 + g1)
            key = (ceil(kappa("dual", k_tfn, a, gm) - EPS), ceil(kappa(model, k_tfn, a, gm) - EPS))
            out[key] = out.get(key, 0.0) + (g1 - g0) * h
    return out


def decision_baseline(key: str, model: str) -> dict:
    bench = BENCHMARKS[INSTANCES[key]]
    lib = library(key)
    levels = range(min(lib) - 1, max(lib) + 1)
    areas = joint_areas(bench["K_tfn"], model, levels)
    assert abs(sum(areas.values()) - 1.0) < 1e-9
    fail = over = same = no_target = 0.0
    excess = []
    for (t, s), a in areas.items():
        assert t in lib and s in lib, (t, s)
        if not lib[t]["feasible"]:
            no_target += a
            continue
        if not lib[s]["feasible"] or not meets(bench, lib[s]["chromosome"], t):
            fail += a
            continue
        e = (lib[s]["cost"] - lib[t]["cost"]) / lib[t]["cost"]
        assert e >= -1e-12, "a design meeting the target cannot be cheaper than the optimum"
        if e > 1e-12:
            over += a
            excess.append((e, a))
        else:
            same += a
    feasible = 1.0 - no_target
    return {
        "fails_target": fail / feasible,
        "meets_at_higher_cost": over / feasible,
        "optimal_cost": same / feasible,
        "max_cost_excess": max((e for e, _ in excess), default=0.0),
        "mean_cost_excess_where_met": sum(e * a for e, a in excess) / (over + same) if over + same else 0.0,
        "area_without_feasible_design": no_target,
    }


# ----------------------------------------------------------------------------- probabilistic counterpart
def cdf_tri(x, a, b, c):
    if x <= a:
        return 0.0
    if x >= c:
        return 1.0
    if x <= b:
        return (x - a) ** 2 / ((c - a) * (b - a))
    return 1.0 - (c - x) ** 2 / ((c - a) * (c - b))


def pdf_tri(x, a, b, c):
    if x < a or x > c:
        return 0.0
    return 2 * (x - a) / ((c - a) * (b - a)) if x <= b else 2 * (c - x) / ((c - a) * (c - b))


def threshold_pmf(k_tfn, dist: str) -> dict[int, float]:
    """q_T = P(ceil(K / M) = T) for K, M independent on the supports of the fuzzy data."""
    k1, k2, k3 = k_tfn
    m1, m2, m3 = 1.0 - SPREAD, 1.0, 1.0 + SPREAD
    if dist == "triangular":
        FK = lambda x: cdf_tri(x, k1, k2, k3)               # noqa: E731
        fM = lambda m: pdf_tri(m, m1, m2, m3)               # noqa: E731
        breaks = [m2]
    elif dist == "uniform":
        FK = lambda x: min(1.0, max(0.0, (x - k1) / (k3 - k1)))   # noqa: E731
        fM = lambda m: 1.0 / (m3 - m1)                            # noqa: E731
        breaks = []
    else:
        raise ValueError(dist)

    def F(L):  # P(K / M <= L) = E[F_K(L M)]
        pts = sorted(set(breaks + [k / L for k in (k1, k2, k3) if m1 < k / L < m3]))
        val, _ = integrate.quad(lambda m: FK(L * m) * fM(m), m1, m3, points=pts or None,
                                epsabs=1e-13, epsrel=1e-12, limit=200)
        return val

    lo, hi = floor(k1 / m3), ceil(k3 / m1)
    cdf = {L: F(L) for L in range(lo, hi + 1)}
    pmf = {T: cdf[T] - cdf[T - 1] for T in range(lo + 1, hi + 1) if cdf[T] - cdf[T - 1] > 0}
    assert abs(sum(pmf.values()) - 1.0) < 1e-9
    return pmf


def tail_matrix(bench: dict, t_max: int):
    """For every design: cost and R(x, T) = P(sum w_i X_i >= T) for T = 0..t_max (exact convolution)."""
    choices = bench["choices"]
    designs = list(itertools.product(*[range(len(c)) for c in choices]))
    costs = np.empty(len(designs), dtype=np.int64)
    tails = np.empty((len(designs), t_max + 1))
    for n, d in enumerate(designs):
        pmf = np.zeros(t_max + 1)
        pmf[0] = 1.0
        for i, j in enumerate(d):
            p, c, w = choices[i][j]
            if j == 0 or w == 0:
                continue
            w = int(w)
            nxt = pmf * (1.0 - p)
            shifted = np.zeros_like(pmf)
            shifted[w:] = pmf[:t_max + 1 - w] if w <= t_max else 0.0
            shifted[t_max] += pmf[t_max + 1 - w:].sum() if w <= t_max else pmf.sum()
            nxt += p * shifted
            pmf = nxt
        tails[n] = np.cumsum(pmf[::-1])[::-1]
        costs[n] = cost(bench, d)
    return designs, costs, tails


def probabilistic(key: str, dist: str, cache: dict) -> dict:
    bench = BENCHMARKS[INSTANCES[key]]
    lib = library(key)
    pmf = threshold_pmf(bench["K_tfn"], dist)
    t_max = max(max(pmf), max(lib))
    if key not in cache:
        cache[key] = tail_matrix(bench, t_max)
    designs, costs, tails = cache[key]
    Ts = np.array(sorted(pmf))
    q = np.array([pmf[t] for t in Ts])
    r_prob = tails[:, Ts] @ q
    ok = r_prob + EPS >= bench["R_target"]
    z = costs[ok].min()
    cand = np.where(ok & (costs == z))[0]
    best = cand[np.argmax(r_prob[cand])]
    design = list(designs[best])
    # spot check of the enumeration against the reliability kernel used everywhere else
    for t in (int(Ts[0]), int(Ts[len(Ts) // 2]), int(Ts[-1])):
        assert abs(tails[best, t] - static_reliability(units(bench, design), t)) < 1e-12
    # the parametric optimum at every T is reproduced by the enumeration
    for t, r in lib.items():
        if r["feasible"]:
            feas = tails[:, t] + EPS >= bench["R_target"]
            assert costs[feas].min() == r["cost"], (key, t)
    t_equiv = max(t for t in lib if meets(bench, design, t))
    in_library = any(r["feasible"] and r["chromosome"] == design for r in lib.values())
    return {"distribution": dist, "design": design, "cost": int(z), "R_prob": float(r_prob[best]),
            "largest_threshold_met": t_equiv, "in_parametric_library": in_library,
            "threshold_support": [int(Ts[0]), int(Ts[-1])],
            "threshold_mean": float(q @ Ts)}


def fixed_design_on_square(key: str, design) -> dict:
    """Fail / higher-cost shares of one fixed design evaluated at every posture of the square."""
    bench = BENCHMARKS[INSTANCES[key]]
    lib = library(key)
    areas = joint_areas(bench["K_tfn"], "crisp", range(min(lib) - 1, max(lib) + 1))  # second index constant
    fail = over = same = no_target = 0.0
    excess = 0.0
    c = cost(bench, design)
    for (t, _), a in areas.items():
        if not lib[t]["feasible"]:
            no_target += a
        elif not meets(bench, design, t):
            fail += a
        elif c > lib[t]["cost"]:
            over += a
            excess = max(excess, (c - lib[t]["cost"]) / lib[t]["cost"])
        else:
            same += a
    feasible = 1.0 - no_target
    return {"fails_target": fail / feasible, "meets_at_higher_cost": over / feasible,
            "optimal_cost": same / feasible, "max_cost_excess": excess}


def rts_necessity(res: dict) -> dict:
    """R_N, R_Pi and crisp adequacy of the RTS-79 probabilistic design, and the cost--necessity frontier
    (exhaustive, rts_true_frontier_v9.json) at the same budget."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("menu", ROOT / "experiments" / "policy_menu_v8.py")
    menu = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(menu)
    r_n, r_p = menu.measures(res["design"])
    bench = BENCHMARKS["power_rts_adapted"]
    r_modal = static_reliability(units(bench, res["design"]), ceil(bench["K_tfn"][1]))
    front = json.loads((RES / "rts_true_frontier_v9.json").read_text())
    best = lambda rows: max((r["r_n"] for r in rows if r["cost"] <= res["cost"]), default=0.0)  # noqa: E731
    return {"R_N": r_n, "R_Pi": r_p, "R_modal": r_modal,
            "frontier_R_N_at_budget": best(front["pareto_frontier"]),
            "adequate_frontier_R_N_at_budget": best(front["adequate_frontier"])}


# ----------------------------------------------------------------------------- driver and table
def main() -> dict:
    out: dict = {"spread": SPREAD, "decision": {}, "probabilistic": {}}
    for key in INSTANCES:
        out["decision"][key] = {m: decision_baseline(key, m) for m in ("crisp", "w_only", "k_only")}
    cache: dict = {}
    for key in PROBABILISTIC:
        out["probabilistic"][key] = {}
        for dist in ("triangular", "uniform"):
            res = probabilistic(key, dist, cache)
            res["on_square"] = fixed_design_on_square(key, res["design"])
            pmf = threshold_pmf(BENCHMARKS[INSTANCES[key]]["K_tfn"], dist)
            res["prob_threshold_above_met"] = sum(v for t, v in pmf.items() if t > res["largest_threshold_met"])
            if key == "rts":
                res["necessity"] = rts_necessity(res)
            out["probabilistic"][key][dist] = res
    (RES / "baselines_v11.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    write_table(out)
    write_supp_tables(out)
    return out


def pct(v: float) -> str:
    return f"{100 * v:.1f}"


def write_table(out: dict) -> None:
    names = {"crisp": "Crisp (modal data)", "w_only": r"Fuzzy capacities only",
             "k_only": r"Fuzzy demand only"}
    cols = ("rts", "fyffe", "pharma")
    lines = [
        r"\begin{table}[!htbp]",
        r"\centering",
        r"\caption{Model baselines evaluated under the declared dual-fuzzy rule. Each baseline uses the "
        r"same data and the same posture $(\alpha,\gamma)$ but a simpler model; its design is evaluated at "
        r"the effective threshold of the dual-fuzzy model. Shares are of the part of the policy square on "
        r"which some design meets the reliability target. ``Fails'': the baseline design does not meet the "
        r"target; ``Costlier'': it meets the target at a cost above the optimum; ``Max.\ excess'': largest "
        r"relative cost excess where it meets the target. Probabilistic counterpart: triangular "
        r"distributions on the supports of the fuzzy data (Section~\ref{sec:baselines_def}).}",
        r"\label{tab:baselines}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{l rrr rrr rrr}",
        r"\hline",
        r" & \multicolumn{3}{c}{IEEE RTS-79} & \multicolumn{3}{c}{Fyffe-derived} & "
        r"\multicolumn{3}{c}{Pharmaceutical-derived}\\",
        r"Baseline & Fails & Costlier & \shortstack{Max.\\excess} & Fails & Costlier & \shortstack{Max.\\excess}"
        r" & Fails & Costlier & \shortstack{Max.\\excess}\\",
        r" & (\%) & (\%) & (\%) & (\%) & (\%) & (\%) & (\%) & (\%) & (\%)\\",
        r"\hline",
    ]
    for m in ("crisp", "w_only", "k_only"):
        cells = []
        for k in cols:
            r = out["decision"][k][m]
            cells += [pct(r["fails_target"]), pct(r["meets_at_higher_cost"]), pct(r["max_cost_excess"])]
        lines.append(names[m] + " & " + " & ".join(cells) + r"\\")
    cells = []
    for k in cols:
        if k in out["probabilistic"]:
            r = out["probabilistic"][k]["triangular"]
            s = r["on_square"]
            cells += [pct(s["fails_target"]), pct(s["meets_at_higher_cost"]), pct(s["max_cost_excess"])]
        else:
            cells += ["--", "--", "--"]
    lines.append("Probabilistic counterpart & " + " & ".join(cells) + r"\\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    (ROOT / "manuscript" / "table_baselines.tex").write_text("\n".join(lines), encoding="utf-8")


def fleet(key: str, design) -> str:
    if key != "rts":
        return "(" + ", ".join(map(str, design)) + ")"
    names = {12: "12~MW oil", 20: "20~MW CT", 50: "50~MW hydro", 76: "76~MW coal", 100: "100~MW oil",
             155: "155~MW coal"}
    bench = BENCHMARKS[INSTANCES[key]]
    count: dict[int, int] = {}
    for i, j in enumerate(design):
        if j:
            w = int(bench["choices"][i][j][2])
            count[w] = count.get(w, 0) + 1
    return " + ".join(f"{c}$\\times${names[w]}" for w, c in sorted(count.items(), key=lambda kv: -kv[0]))


def write_supp_tables(out: dict) -> None:
    label = {"rts": "RTS-79", "fyffe": "Fyffe-derived", "pharma": "Pharma-derived"}
    names = {"crisp": "Crisp", "w_only": "Fuzzy capacities only", "k_only": "Fuzzy demand only"}
    a = [r"\begin{table}[!htbp]", r"\centering",
         r"\caption{Decision-level model baselines in full (shares of the part of the policy square on which some "
         r"design meets the target). ``Optimal'': the baseline design is also optimal under the declared rule; "
         r"``Mean excess'': area-weighted relative cost excess over the part on which the baseline design meets the "
         r"target, including the part on which it is optimal.}",
         r"\label{tab:supp_baselines}", r"\footnotesize", r"\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{llrrrrr}", r"\hline",
         r"Instance & Baseline & Fails & Costlier & Optimal & Max.\ excess & Mean excess\\", r" & & (\%) & (\%) & (\%) & (\%) & (\%)\\",
         r"\hline"]
    for k in ("rts", "fyffe", "pharma"):
        for m in ("crisp", "w_only", "k_only"):
            r = out["decision"][k][m]
            a.append(f"{label[k] if m == 'crisp' else ''} & {names[m]} & {pct(r['fails_target'])} & "
                     f"{pct(r['meets_at_higher_cost'])} & {pct(r['optimal_cost'])} & {pct(r['max_cost_excess'])} & "
                     f"{pct(r['mean_cost_excess_where_met'])}\\\\")
        if k != "pharma":
            a.append(r"\hline")
    a += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    b = [r"\begin{table}[!htbp]", r"\centering",
         r"\caption{Probabilistic counterpart: optimal design under triangular and uniform densities on the supports "
         r"of the fuzzy data, the largest effective threshold at which it meets the target, the probability that the "
         r"random threshold exceeds it, and its evaluation over the policy square under the declared rule. Pharmaceutical "
         r"designs are given as option vectors.}",
         r"\label{tab:supp_probabilistic}", r"\footnotesize", r"\setlength{\tabcolsep}{2.5pt}",
         r"\begin{tabular}{llp{3.3cm}rrrrrrr}", r"\hline",
         r"Instance & Density & Design & Cost & $R_{\mathrm{prob}}$ & $T$ met & $\Pr\{\hat{T}>\cdot\}$ & Fails & "
         r"Costlier & Max.\ exc.\\", r" & & & & & & & (\%) & (\%) & (\%)\\", r"\hline"]
    for k in PROBABILISTIC:
        for d in ("triangular", "uniform"):
            r = out["probabilistic"][k][d]
            s_ = r["on_square"]
            b.append(f"{label[k] if d == 'triangular' else ''} & {d} & {fleet(k, r['design'])} & {r['cost']} & "
                     f"{r['R_prob']:.4f} & {r['largest_threshold_met']} & {r['prob_threshold_above_met']:.3f} & "
                     f"{pct(s_['fails_target'])} & {pct(s_['meets_at_higher_cost'])} & {pct(s_['max_cost_excess'])}\\\\")
    b += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    (ROOT / "manuscript" / "table_baselines_supp.tex").write_text("\n".join(a + b), encoding="utf-8")


if __name__ == "__main__":
    res = main()
    for k, v in res["decision"].items():
        for m, r in v.items():
            print(k, m, {a: round(b, 4) for a, b in r.items()})
    for k, v in res["probabilistic"].items():
        for d, r in v.items():
            print(k, d, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in r.items() if a not in ("on_square", "necessity")},
                  {a: round(b, 4) for a, b in r["on_square"].items()}, r.get("necessity"))
