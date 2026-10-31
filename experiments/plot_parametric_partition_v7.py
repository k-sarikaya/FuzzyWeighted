"""Complete parametric solution of the IEEE RTS-79 fleet over the policy square.

Panel (a): the (alpha, gamma) square partitioned by the optimal design.  Thin lines are
           the level curves kappa(alpha, gamma) = T of every attainable integer threshold,
           i.e. the boundaries between the 73 crisp subproblems; thick lines are the
           closed-form boundaries of Corollary 2 at which the optimal design changes;
           markers are the 30 cells of the grid campaign.
Panel (b): the same solution as a function of the single effective threshold kappa:
           minimum cost Z* on each interval (T-1, T], coloured by design, with the kappa
           values reached by the 30 grid cells marked on the axis.

Inputs : results/generated/parametric_stability_rts.json  (T -> certified crisp optimum)
         results/generated/benchmark_power_rts.json       (30-cell grid campaign)
Output : figures/fig_parametric_partition.pdf / .png

The script checks, before plotting, that (i) the closed-form boundary of Corollary 2
lies on the level set kappa = T, and (ii) every grid cell falls in the region whose design
and cost the grid campaign certified.  It stops if either check fails.

    python experiments/plot_parametric_partition_v7.py
"""

from __future__ import annotations

import json
from math import ceil
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "generated"
OUT = ROOT / "figures"

K_TFN = (162.0, 180.0, 198.0)   # demand TFN of the RTS-79 case study (MW)
SPREAD = 0.10                    # common proportional capacity spread

INK, INK2 = "#1f1f1e", "#5f5e59"
# Ordered light -> dark ramp: designs are indexed by increasing cost.
RAMP = ["#fde6a8", "#f5bd72", "#e58e63", "#c46679", "#8b5a9c", "#4b5ca6", "#1e3a72"]

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "savefig.dpi": 300,
})


def kappa(alpha, gamma):
    """Effective threshold of Theorem 1 (vectorised)."""
    k1, k2, k3 = K_TFN
    lower = k1 + alpha * (k2 - k1)
    upper = k3 - alpha * (k3 - k2)
    demand = upper - gamma * (upper - lower)
    return demand / (1.0 + SPREAD * (1.0 - alpha) * (2.0 * gamma - 1.0))


def gamma_boundary(alpha, level):
    """Corollary 2: the gamma at which kappa(alpha, gamma) = level (closed form)."""
    k1, k2, k3 = K_TFN
    lower = k1 + alpha * (k2 - k1)
    upper = k3 - alpha * (k3 - k2)
    with np.errstate(divide="ignore", invalid="ignore"):   # alpha = 1: no boundary
        return (upper - level * (1.0 - SPREAD * (1.0 - alpha))) / (
            (upper - lower) + 2.0 * SPREAD * level * (1.0 - alpha))


def load():
    param = json.loads((RESULTS / "parametric_stability_rts.json").read_text(encoding="utf-8"))
    grid = json.loads((RESULTS / "benchmark_power_rts.json").read_text(encoding="utf-8"))
    regions = [r for r in param["stability_regions"] if r["feasible"]]
    assert len(regions) == len(param["stability_regions"]), "infeasible threshold in the path"
    assert all(regions[i]["T_max"] + 1 == regions[i + 1]["T_min"] for i in range(len(regions) - 1))
    return param, grid, regions


def check(param, grid, regions) -> list[int]:
    """Closed-form boundaries lie on kappa = T; each grid cell is in its certified region."""
    alphas = np.linspace(0.0, 1.0, 401)
    for r in regions[1:]:
        level = r["T_min"] - 1
        g = gamma_boundary(alphas, level)
        inside = (g >= 0) & (g <= 1)
        assert np.allclose(kappa(alphas[inside], g[inside]), level, atol=1e-9)
    hit = []
    for cell in grid:
        T = ceil(kappa(cell["alpha"], cell["gamma"]) - 1e-12)
        idx = next(i for i, r in enumerate(regions) if r["T_min"] <= T <= r["T_max"])
        r = regions[idx]
        assert T == cell["T"], (cell["alpha"], cell["gamma"], T, cell["T"])
        assert cell["design"] == r["chromosome"] and cell["cost"] == r["cost"], cell
        hit.append(idx)
    return hit


def label_point(mask, A, G):
    """Point of the region farthest from its complement (a readable place for a label)."""
    step = max(1, mask.shape[0] // 90)
    m, a, g = mask[::step, ::step], A[::step, ::step], G[::step, ::step]
    inside = np.column_stack([a[m], g[m]])
    outside = np.column_stack([a[~m], g[~m]])
    best, best_d = inside[0], -1.0
    for chunk in np.array_split(inside, max(1, len(inside) // 500)):
        d = np.sqrt(((chunk[:, None, :] - outside[None, :, :]) ** 2).sum(-1)).min(1)
        # the square's edges also bound the region
        d = np.minimum(d, np.minimum.reduce([chunk[:, 0], 1 - chunk[:, 0], chunk[:, 1], 1 - chunk[:, 1]]) + 0.02)
        j = int(d.argmax())
        if d[j] > best_d:
            best, best_d = chunk[j], d[j]
    return best, best_d


def main() -> None:
    param, grid, regions = load()
    hit = check(param, grid, regions)
    n = len(regions)
    k_min, k_max = param["kappa_min"], param["kappa_max"]
    T_all = param["integer_thresholds"]

    N = 801
    A, G = np.meshgrid(np.linspace(0, 1, N), np.linspace(0, 1, N))
    K = kappa(A, G)
    # region i holds kappa in (T_min_i - 1, T_max_i]
    edges = [k_min - 1e-9] + [r["T_max"] for r in regions]
    D = np.digitize(K, edges[1:-1], right=True)
    area = [float((D == i).mean()) for i in range(n)]

    cmap = ListedColormap(RAMP[:n])
    norm = BoundaryNorm(np.arange(-0.5, n), n)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.35),
                                   gridspec_kw={"width_ratios": [1.0, 1.25], "wspace": 0.28})

    # ---- (a) partition of the policy square ------------------------------------------
    ax1.pcolormesh(A, G, D, cmap=cmap, norm=norm, shading="auto", rasterized=True)
    ax1.contour(A, G, K, levels=[t for t in T_all if k_min < t < k_max],
                colors="white", linewidths=0.25, alpha=0.55)
    alpha_fine = np.linspace(0, 1, 1200)
    for r in regions[1:]:
        g = gamma_boundary(alpha_fine, r["T_min"] - 1)
        ok = (g >= 0) & (g <= 1)
        ax1.plot(alpha_fine[ok], g[ok], color=INK, linewidth=1.1)
    ga = np.array([c["alpha"] for c in grid])
    gg = np.array([c["gamma"] for c in grid])
    ax1.scatter(ga, gg, s=14, facecolor="white", edgecolor=INK, linewidth=0.7, zorder=5,
                clip_on=False, label="30 grid cells")
    for i in range(n):
        (pa, pg), depth = label_point(D == i, A, G)
        dark = i >= 4
        if depth > 0.035:
            ax1.text(pa, pg, f"D{i + 1}", ha="center", va="center", fontsize=8, fontweight="bold",
                     color="white" if dark else INK, zorder=6)
        else:   # region too thin for an inside label: label above the square, leader line
            ax1.annotate(f"D{i + 1}", (pa, pg), xytext=(0.07 + 0.16 * i, 1.035),
                         textcoords="data", fontsize=8, fontweight="bold", color=INK,
                         ha="center", va="bottom", annotation_clip=False, zorder=7,
                         arrowprops=dict(arrowstyle="-", color=INK, linewidth=0.6,
                                         shrinkA=0, shrinkB=0))
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_aspect("equal")
    ax1.set_xlabel(r"cut level $\alpha$")
    ax1.set_ylabel(r"optimism index $\gamma$")
    ax1.set_title(r"(a) Policy square $(\alpha,\gamma)$", loc="left", pad=16)
    ax1.legend(loc="lower right", frameon=True, framealpha=0.9, handletextpad=0.3,
               borderpad=0.3, fontsize=7)

    # ---- (b) the same solution as a function of kappa --------------------------------
    for i, r in enumerate(regions):
        lo = max(r["T_min"] - 1, k_min)
        hi = r["T_max"]
        ax2.axvspan(lo, hi, color=RAMP[i], alpha=0.45, linewidth=0)
        ax2.hlines(r["cost"], lo, hi, color=INK, linewidth=1.6)
        note = "" if i in hit else " (no grid cell)"
        text = f"D{i + 1}  \\${r['cost']}M{note}"
        if hi - lo < 8:     # narrow band: vertical label inside the band
            above = r["cost"] < 170          # tall labels of high steps go below the step
            ax2.text((lo + hi) / 2, r["cost"] + (4 if above else -4), text, rotation=90,
                     ha="center", va="bottom" if above else "top", fontsize=7, color=INK)
        else:
            ax2.text((lo + hi) / 2, r["cost"] + 3, text, ha="center", va="bottom",
                     fontsize=7, color=INK)
        if i:
            ax2.vlines(lo, regions[i - 1]["cost"], r["cost"], color=INK, linewidth=0.8,
                       linestyles=":")
    gk = np.array([kappa(c["alpha"], c["gamma"]) for c in grid])
    ax2.plot(gk, np.full_like(gk, 122), linestyle="none", marker="|", markersize=9,
             color=INK, markeredgewidth=0.9, label=r"$\kappa$ of the 30 grid cells")
    ax2.set_xlim(k_min, k_max)
    ax2.set_ylim(118, 222)
    ax2.set_xlabel(r"effective threshold $\kappa(\alpha,\gamma)$ [MW]")
    ax2.set_ylabel(r"minimum cost $Z^*$ [M\$]")
    ax2.set_title(f"(b) {len(T_all)} crisp problems, {n} optimal designs", loc="left", pad=20)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(loc="upper left", frameon=False, fontsize=7)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig_parametric_partition.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_parametric_partition.png", bbox_inches="tight")
    plt.close(fig)

    missed = [i + 1 for i in range(n) if i not in hit]
    print(f"kappa in [{k_min:.2f}, {k_max:.2f}]; {len(T_all)} thresholds; {n} designs")
    for i, r in enumerate(regions):
        print(f"  D{i + 1}: T {r['T_min']}-{r['T_max']}  cost {r['cost']}  "
              f"R {r['reliability']:.4f}  area {100 * area[i]:.1f}%  "
              f"{'on grid' if i in hit else 'NOT on grid'}")
    print(f"grid cells: {len(grid)} checked, designs reached {sorted(set(h + 1 for h in hit))}, "
          f"missed {missed}")


if __name__ == "__main__":
    main()
