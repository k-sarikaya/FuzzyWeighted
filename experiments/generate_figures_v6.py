"""Every figure of manuscript_v6, rebuilt from the v6 artifacts.

Inputs (no solver is run here):
  results/generated/benchmark_tables_v6.json   exact optimum of every grid cell
  results/generated/benchmark_power_rts.json   exact optimum of every IEEE RTS-79 cell

Outputs in figures/ (PDF for the manuscript, PNG for quick viewing):
  Main text
  fig2_cost_sensitivity  RTS-79 minimum cost against alpha, one line per gamma
  fig5_policy_effect     RTS-79 cost and necessity reliability against gamma
  Supplementary Material (Fyffe benchmark)
  figS_fyffe_sensitivity Fyffe cost against alpha; cost and necessity reliability against gamma
  fig3_cost_heatmap      Fyffe minimum cost over the grid
  fig4_config_map        Fyffe number of installed subsystems over the grid

    python3 experiments/generate_figures_v6.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "generated"
OUT = ROOT / "figures"

INK, INK2, GRID = "#1f1f1e", "#5f5e59", "#e4e3de"
# Diverging pessimistic (red) -> neutral (gray) -> optimistic (blue) for gamma.
GAMMA_STYLE = {
    0.0: ("#c0392f", "-", "o"),
    0.3: ("#e8908b", "-", "D"),
    0.5: ("#6f6e69", "--", "s"),
    0.7: ("#6da7ec", "-", "v"),
    1.0: ("#1c5cab", "-", "^"),
}
# Ordinal blue ramp for alpha (validated: monotone L, visible gaps).
ALPHA_RAMP = {0.0: "#86b6ef", 0.3: "#3987e5", 0.5: "#256abf", 0.7: "#184f95", 1.0: "#0d366b"}
# Three categorical slots for the three evaluation cells of Figure 1.
CELL_STYLE = [("#eb6834", "-", "o"), ("#6f6e69", "--", "s"), ("#2a78d6", "-.", "^")]
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
NAMES = {"pharma_adapted": "Pharmaceutical Plant (adapted)", "fyffe_adapted": "Fyffe 14-subsystem (adapted)",
         "power_rts_adapted": "IEEE RTS-79 fleet (cost in M$)"}
SHORT = {"pharma_adapted": "Pharmaceutical Plant", "fyffe_adapted": "Fyffe 14-subsystem"}
SHORT_ALL = {**SHORT, "power_rts_adapted": "IEEE RTS-79"}
# Label offsets (points) for the T* annotations, chosen to avoid the curves.
T_STAR_OFFSET = {"pharma_adapted": [(6, -12), (6, -12), (6, -12)],
                 "fyffe_adapted": [(6, -12), (-8, 20), (6, -12)]}

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300,
})


def _rows(tables: dict, catalog: str) -> list[dict]:
    block = tables[catalog]
    return block["rows"] if isinstance(block, dict) else block


def _cell(rows: list[dict], alpha: float, gamma: float) -> dict:
    return next(r for r in rows if abs(r["alpha"] - alpha) < 1e-9 and abs(r["gamma"] - gamma) < 1e-9)


def _grid(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def _save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def _cost_vs_alpha(ax, tables: dict, key: str, legend: str = "upper right") -> None:
    rows = _rows(tables, key)
    if True:
        alphas = sorted({r["alpha"] for r in rows})
        for gamma, (color, style, marker) in GAMMA_STYLE.items():
            values = [(_cell(rows, a, gamma)["cost"] if _cell(rows, a, gamma)["feasible"] else np.nan)
                      for a in alphas]
            ax.plot(alphas, values, color=color, linestyle=style, marker=marker, markersize=5,
                    linewidth=1.6, label=rf"$\gamma={gamma:g}$")
            for index, (a, v) in enumerate(zip(alphas, values)):
                if np.isnan(v):
                    anchor = next(x for x in values[index:] if not np.isnan(x))
                    ax.plot([a], [anchor], marker="x", color=color, markersize=6)
                    ax.annotate("infeasible", (a, anchor), xytext=(6, -3),
                                textcoords="offset points", fontsize=7, color=INK2)
        ax.set_title(NAMES[key])
        ax.set_xlabel(r"$\alpha$")
        ax.set_ylabel("minimum feasible cost (exact)")
        ax.set_xticks(alphas)
        _grid(ax)
        if legend == "outside":
            ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
        else:
            ax.legend(frameon=False, loc=legend)


def fig2(tables: dict) -> None:
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    _cost_vs_alpha(ax, tables, "power_rts_adapted", legend="outside")
    fig.tight_layout()
    _save(fig, "fig2_cost_sensitivity")


def _heatmap(tables: dict, key: str, value, label: str, name: str) -> None:
    rows = _rows(tables, key)
    alphas = sorted({r["alpha"] for r in rows})
    gammas = sorted({r["gamma"] for r in rows})
    grid = np.array([[value(_cell(rows, a, g)) for a in alphas] for g in gammas], dtype=float)
    fig, ax = plt.subplots(figsize=(3.5, 2.9))
    image = ax.imshow(grid, origin="lower", cmap=SEQ, aspect="auto")
    lo, hi = np.nanmin(grid), np.nanmax(grid)
    for i in range(len(gammas)):
        for j in range(len(alphas)):
            v = grid[i, j]
            dark = (v - lo) / (hi - lo if hi > lo else 1) > 0.55
            ax.text(j, i, f"{v:g}", ha="center", va="center", fontsize=8,
                    color="white" if dark else INK)
    ax.set_xticks(range(len(alphas)), [f"{a:g}" for a in alphas])
    ax.set_yticks(range(len(gammas)), [f"{g:g}" for g in gammas])
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\gamma$")
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.colorbar(image, ax=ax, fraction=0.05, pad=0.03)
    bar.set_label(label)
    bar.outline.set_visible(False)
    fig.tight_layout()
    _save(fig, name)


def _policy_effect(row_axes, tables: dict, key: str) -> None:
    rows = _rows(tables, key)
    if True:
        gammas = sorted({r["gamma"] for r in rows})
        for alpha, color in ALPHA_RAMP.items():
            cells = [_cell(rows, alpha, g) for g in gammas]
            row_axes[0].plot(gammas, [c["cost"] for c in cells], color=color, marker="o", markersize=4,
                             linewidth=1.6, label=rf"$\alpha={alpha:g}$")
            row_axes[1].plot(gammas, [c["R_necessity"] for c in cells], color=color, marker="o", markersize=4,
                             linewidth=1.6, label=rf"$\alpha={alpha:g}$")
        row_axes[0].set_ylabel("minimum feasible cost (exact)")
        row_axes[1].set_ylabel(r"necessity reliability $R_N$")
        row_axes[1].set_ylim(-0.02, 1.0)
        for ax, title in zip(row_axes, ("Cost", "Strength of the guarantee")):
            ax.set_xlabel(r"$\gamma$")
            ax.set_xticks(gammas)
            ax.set_title(f"{SHORT_ALL[key]}: {title.lower()}")
            _grid(ax)
    row_axes[0].legend(frameon=False, loc="lower left", fontsize=7)


def fig5(tables: dict) -> None:
    """Cost and necessity reliability against gamma on the RTS-79 fleet."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
    _policy_effect(axes, tables, "power_rts_adapted")
    fig.tight_layout()
    _save(fig, "fig5_policy_effect")


def fig_s_fyffe(tables: dict) -> None:
    """Supplementary: the Fyffe counterparts of Figures 2 and 5 in one figure."""
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.0))
    _cost_vs_alpha(axes[0], tables, "fyffe_adapted")
    axes[0].set_title("Fyffe 14-subsystem: cost against $\\alpha$")
    _policy_effect(axes[1:], tables, "fyffe_adapted")
    fig.tight_layout()
    _save(fig, "figS_fyffe_sensitivity")


def main() -> None:
    tables = json.loads((RESULTS / "benchmark_tables_v6.json").read_text())
    tables["power_rts_adapted"] = {"rows": json.loads((RESULTS / "benchmark_power_rts.json").read_text())}
    fig2(tables)
    _heatmap(tables, "fyffe_adapted", lambda r: r["cost"] if r["feasible"] else np.nan,
             "minimum feasible cost", "fig3_cost_heatmap")
    _heatmap(tables, "fyffe_adapted", lambda r: sum(1 for c in r["design"] if c != 0) if r["feasible"] else np.nan,
             "installed subsystems", "fig4_config_map")
    fig5(tables)
    fig_s_fyffe(tables)
    print("figures written to", OUT)


if __name__ == "__main__":
    main()
