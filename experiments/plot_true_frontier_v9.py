"""Plot True Cost-R_N Pareto Frontier vs Hurwicz Policy Menu for IEEE RTS-79 fleet.

Resolves Reviewer concerns R2-M3 and R3-M3:
- Visualizes the full catalog (3^10 = 59,049 designs) in the Cost - R_N plane.
- Traces the true Pareto-optimal boundary (min Cost s.t. R_N >= lambda).
- Overlays the 7 Hurwicz-optimal fleet designs (D1..D7) from the policy menu.
- Traces the adequacy-constrained frontier (additionally R_modal >= 0.95 at the modal demand).
- Marks the $171M catalog fleet with the largest R_N (dominated by the $170M frontier design).
"""

import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "generated"
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

frontier_data = json.loads((RESULTS / "rts_true_frontier_v9.json").read_text(encoding="utf-8"))
pareto = frontier_data["pareto_frontier"]
cost_171 = frontier_data["cost_171_designs"]

# Hurwicz menu D1..D7, read from the generated policy menu (no hand-typed values)
policy = json.loads((RESULTS / "policy_menu_v8.json").read_text(encoding="utf-8"))
menu = [{"label": d["label"], "cost": d["cost"], "r_n": d["R_N"]} for d in policy["menu"]]
adequate = [p for p in frontier_data["adequate_frontier"] if p["cost"] <= 240]

# Set publication style
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9,
    "figure.titlesize": 12
})

fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=300)

# True Pareto Frontier coordinates (stepped)
p_costs = [p["cost"] for p in pareto if p["cost"] > 0]
p_rns = [p["r_n"] for p in pareto if p["cost"] > 0]

# Unconstrained frontier (min Z s.t. R_N >= lambda)
ax.step(p_costs, p_rns, where="post", color="#1f77b4", linewidth=2.0,
        label=r"Unconstrained frontier (all $3^{10}$ designs)", zorder=3)
ax.scatter(p_costs, p_rns, color="#1f77b4", s=25, zorder=4)

# Adequacy-constrained frontier (additionally R_modal >= 0.95)
a_costs = [p["cost"] for p in adequate]
a_rns = [p["r_n"] for p in adequate]
ax.step(a_costs, a_rns, where="post", color="#ff7f0e", linewidth=2.0, linestyle="--",
        label=r"Adequacy-constrained frontier ($R_{\mathrm{modal}} \geq 0.95$)", zorder=3)
ax.scatter(a_costs, a_rns, color="#ff7f0e", s=30, marker="s", zorder=4)

# Plot Hurwicz Menu points
m_costs = [m["cost"] for m in menu]
m_rns = [m["r_n"] for m in menu]
ax.scatter(m_costs, m_rns, color="#d62728", s=65, marker="D", edgecolor="black", linewidth=0.8, label="Hurwicz Policy Menu (D1–D7)", zorder=5)

offsets = {"D1": (-3, 0.04), "D2": (1, -0.065), "D3": (-3, 0.04), "D4": (-3, -0.065),
           "D5": (1, -0.07), "D6": (-15, 0.04), "D7": (-15, 0.04)}
for m in menu:
    dx, dy = offsets[m["label"]]
    ax.annotate(f"{m['label']} (${m['cost']}M)", (m["cost"], m["r_n"]),
                xytext=(m["cost"] + dx, m["r_n"] + dy),
                fontsize=8.5, fontweight="bold", color="#d62728")

# $171M catalog fleet with the largest R_N (dominated by the $170M frontier design)
best_171 = max(cost_171, key=lambda x: x["r_n"])
ax.scatter(best_171["cost"], best_171["r_n"], color="#2ca02c", s=70, marker="^", edgecolor="black",
           linewidth=0.8, zorder=6, label=f"\\$171M catalog fleet ($R_N={best_171['r_n']:.3f}$, dominated)")

# Formatting
ax.set_xlabel("Capital Investment Cost [M$]", fontweight="bold")
ax.set_ylabel(r"Necessity Reliability $R_N$ of Meeting Demand", fontweight="bold")
ax.set_xlim(120, 225)
ax.set_ylim(-0.12, 1.05)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(loc="lower right", framealpha=0.95)

plt.tight_layout()
fig_pdf = FIGURES / "fig_true_pareto_frontier.pdf"
fig_png = FIGURES / "fig_true_pareto_frontier.png"
fig.savefig(fig_pdf, bbox_inches="tight")
fig.savefig(fig_png, bbox_inches="tight")
plt.close(fig)

print(f"Generated {fig_pdf} and {fig_png}")
