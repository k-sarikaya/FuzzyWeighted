"""Crisp versus dual-fuzzy reliability of the same designs (IEEE RTS-79 fleet).

For each of the seven designs of the complete parametric solution, the static reliability is
evaluated, for the same installed units, under five success conditions:

  crisp                     modal capacities against the modal demand,   W_S >= k2
  necessity cell (0, 0):    capacities fuzzy only      W_S (1 - d) >= k2
                            threshold fuzzy only       W_S         >= k3
                            both fuzzy (the model)     W_S (1 - d) >= k3
  optimistic boundary cell (0, 1):  capacities fuzzy only      W_S (1 + d) >= k2
                            threshold fuzzy only       W_S         >= k1
                            both fuzzy (the model)     W_S (1 + d) >= k1

with d = 0.10 and k = (162, 180, 198) MW.  At alpha = 0 the Hurwicz rule compares the endpoints of
the support, so the "both" columns are the reliabilities of the model at those cells; the crisp
value is the Hurwicz reliability at alpha = 1 (and at gamma = 1/2 for this symmetric data).
Exact rational arithmetic, enumeration of all up/down states of the installed units.

Writes results/generated/crisp_vs_fuzzy_v8.json and manuscript/table_crisp_vs_fuzzy.tex.

    python experiments/crisp_vs_fuzzy_v8.py
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuzzy_wkon.benchmarks import BENCHMARKS  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "campaign", ROOT / "experiments" / "run_power_rts_campaign_v7.py")
campaign = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(campaign)

RESULTS = ROOT / "results" / "generated"
B = BENCHMARKS["power_rts_adapted"]
K1, K2, K3 = (F(str(v)) for v in B["K_tfn"])
D = F(1, 10)

# (capacity factor, threshold) for each condition
CONDITIONS = {
    "crisp": (F(1), K2),
    "nec_w_only": (1 - D, K2),
    "nec_k_only": (F(1), K3),
    "nec_both": (1 - D, K3),
    "pos_w_only": (1 + D, K2),
    "pos_k_only": (F(1), K1),
    "pos_both": (1 + D, K1),
}


def units(chromosome: list[int]) -> list[tuple[F, F]]:
    """(modal capacity, availability) of every installed unit."""
    out = []
    for choices, c in zip(B["choices"], chromosome):
        if c:
            p, _cost, w = choices[c]
            out.append((F(str(w)), F(str(p))))
    return out


def reliability(design_units: list[tuple[F, F]], factor: F, threshold: F) -> F:
    total = F(0)
    for states in itertools.product((0, 1), repeat=len(design_units)):
        w = sum(u[0] for u, s in zip(design_units, states) if s)
        if factor * w >= threshold:
            p = F(1)
            for (_, q), s in zip(design_units, states):
                p *= q if s else 1 - q
            total += p
    return total


def main() -> None:
    param = json.loads((RESULTS / "parametric_stability_rts.json").read_text(encoding="utf-8"))
    regions = [r for r in param["stability_regions"] if r["feasible"]]
    rows = []
    for i, r in enumerate(regions):
        u = units(r["chromosome"])
        fleet, mw, _ = campaign.describe_fleet(B, r["chromosome"])
        vals = {name: reliability(u, f, k) for name, (f, k) in CONDITIONS.items()}
        rows.append({"label": f"D{i + 1}", "fleet": fleet, "capacity_mw": mw, "cost": r["cost"],
                     "chromosome": r["chromosome"], **{k: float(v) for k, v in vals.items()}})

    out = {"description": __doc__.split("\n\n")[0], "spread": float(D),
           "threshold_tfn": [float(K1), float(K2), float(K3)], "rows": rows}
    (RESULTS / "crisp_vs_fuzzy_v8.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    def f4(x: float) -> str:
        return "0" if x == 0 else f"{x:.4f}"

    body = "\n".join(
        f"{r['label']} & {r['capacity_mw']} & {r['cost']} & {f4(r['crisp'])} & "
        f"{f4(r['nec_w_only'])} & {f4(r['nec_k_only'])} & {f4(r['nec_both'])} & "
        f"{f4(r['pos_w_only'])} & {f4(r['pos_k_only'])} & {f4(r['pos_both'])} \\\\" for r in rows)
    table = (
        "\\begin{table}[!htbp]\n\\providecommand{\\revtablecolor}{}\\revtablecolor\n\\centering\n"
        "\\caption{Crisp and dual-fuzzy reliability of the same designs (IEEE RTS-79 fleet, designs D1--D7 of "
        "Figure~\\ref{fig:parametric}). ``Crisp'': modal capacities against the modal demand of "
        "$180$~MW. In the necessity cell $(0,0)$ and the optimistic boundary cell $(0,1)$, ``$\\tilde{w}$ only'' "
        "fuzzifies the capacities and keeps the demand at $180$~MW, ``$\\tilde{k}$ only'' fuzzifies the "
        "demand and keeps the modal capacities, and ``both'' is the dual-fuzzy model. "
        "}\n"
        "\\label{tab:crisp_vs_fuzzy}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{lrrrrrrrrr}\n\\hline\n"
        " & & & & \\multicolumn{3}{c}{Necessity cell $(0,0)$} & \\multicolumn{3}{c}{Optimistic boundary cell $(0,1)$}\\\\\n"
        "Design & MW & Cost (M\\$) & Crisp & $\\tilde{w}$ only & $\\tilde{k}$ only & both & "
        "$\\tilde{w}$ only & $\\tilde{k}$ only & both\\\\\n\\hline\n"
        + body + "\n\\hline\n\\end{tabular}\n\\end{table}\n")
    (ROOT / "manuscript" / "table_crisp_vs_fuzzy.tex").write_text(table, encoding="utf-8")

    for r in rows:
        print(f"{r['label']} {r['capacity_mw']:4d} MW {r['cost']:4d} crisp={r['crisp']:.4f} | "
              f"nec w={r['nec_w_only']:.4f} k={r['nec_k_only']:.4f} both={r['nec_both']:.4f} | "
              f"pos w={r['pos_w_only']:.4f} k={r['pos_k_only']:.4f} both={r['pos_both']:.4f}")


if __name__ == "__main__":
    main()
