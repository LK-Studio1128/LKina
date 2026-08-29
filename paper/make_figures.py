#!/usr/bin/env python3
"""Generate the three figures for the LKina manuscript (v2: cleaner layout)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/figures"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.edgecolor": "#444",
    "axes.linewidth": 0.9,
    "axes.titlesize": 12,
    "axes.labelsize": 11.5,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 200,
})

C_AD4   = "#9AA5B1"
C_AD4ZN = "#E8A33D"
C_LKINA = "#C0392B"
C_VINA  = "#9AA5B1"
C_P1    = "#E8A33D"
C_P1P2  = "#C0392B"

# ----------------------------------------------------------------------------
# Figure 1 — Metal redocking success rates
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.8, 4.6))

groups = ["Zn²⁺\n(292 complexes)", "Fe³⁺\n(56 complexes)", "Cu²⁺\n(41 complexes)"]
methods = {
    "Standard AD4":       [58.2, 48.2, 39.0],
    "AutoDock4Zn":        [70.5, None, None],
    "LKina (TZ/MH/JT)":   [74.3, 66.1, 63.4],
}
colors = {"Standard AD4": C_AD4, "AutoDock4Zn": C_AD4ZN, "LKina (TZ/MH/JT)": C_LKINA}

x = np.arange(len(groups))
n = len(methods)
width = 0.25
offsets = np.linspace(-(n - 1) / 2, (n - 1) / 2, n) * width

for i, (name, vals) in enumerate(methods.items()):
    xs = x + offsets[i]
    ys = [v if v is not None else 0 for v in vals]
    bars = ax.bar(xs, ys, width, label=name, color=colors[name],
                  edgecolor="white", linewidth=0.6, zorder=3)
    for xi, yi, v in zip(xs, ys, vals):
        if v is not None:
            ax.text(xi, yi + 1.4, f"{v:.1f}", ha="center", va="bottom",
                    fontsize=9.5, fontweight="bold", color=colors[name])

ax.set_xticks(x)
ax.set_xticklabels(groups)
ax.set_ylabel("Top-1 RMSD ≤ 2.0 Å success rate (%)")
ax.set_ylim(0, 86)
ax.axhline(0, color="#444", linewidth=0.9)
ax.grid(axis="y", color="#E5E7EB", linewidth=0.8, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1.0, 1.12))

# Improvements over AD4 — placed as a small caption box in the upper-left
ax.text(0.02, 1.12,
        "Improvement over standard AD4:\n"
        "  Zn²⁺   +16.1 pp\n"
        "  Fe³⁺   +17.9 pp\n"
        "  Cu²⁺   +24.4 pp",
        transform=ax.transAxes, fontsize=9.5, va="bottom", ha="left",
        color=C_LKINA, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#FCEBEB",
                  edgecolor=C_LKINA, linewidth=0.8))

fig.tight_layout()
fig.savefig(f"{OUT}/fig1_metal_success.png", bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------------------------
# Figure 2 — Cumulative top-1 RMSD distribution (Zn benchmark)
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.8, 4.4))

rmsd = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
cum = {
    "Standard AD4": [16.1, 31.5, 53.8, 58.2, 64.7, 70.2, 73.6, 76.4],
    "AutoDock4Zn":  [26.0, 43.8, 65.4, 70.5, 76.4, 80.8, 83.6, 85.6],
    "LKina TZ":     [29.5, 46.9, 67.1, 74.3, 80.1, 84.6, 87.3, 89.0],
}
styles = {"Standard AD4": ("--", C_AD4), "AutoDock4Zn": ("-.", C_AD4ZN), "LKina TZ": ("-", C_LKINA)}

for name, (ls, c) in styles.items():
    ax.plot(rmsd, cum[name], ls, color=c, linewidth=2.2, marker="o",
            markersize=5.5, label=name, zorder=3)

ax.axvline(2.0, color="#666", linestyle=":", linewidth=1.0)
ax.text(2.0, 12, "2.0 Å", color="#666", fontsize=9.5, ha="center")
ax.set_xlabel("Top-1 RMSD threshold (Å)")
ax.set_ylabel("Cumulative success rate (%)")
ax.set_ylim(0, 100)
ax.set_xlim(0.4, 4.1)
ax.grid(axis="y", color="#E5E7EB", linewidth=0.8, zorder=0)
ax.grid(axis="x", color="#EEF0F3", linewidth=0.6, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, loc="lower right")

fig.tight_layout()
fig.savefig(f"{OUT}/fig2_zn_cumulative_rmsd.png", bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------------------------
# Figure 3 — Covalent docking benchmark
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.8, 4.2))

labels = ["Standard\nVina", "LKina\nP1 (distance)", "LKina\nP1+P2\n(distance+angle)"]
rmsd_ok   = [42.9, 60.7, 64.3]
nac_ok    = [None, 67.9, 71.4]
colors3   = [C_VINA, C_P1, C_P1P2]

x = np.arange(len(labels))
bars = ax.bar(x, rmsd_ok, 0.55, color=colors3, edgecolor="white",
              linewidth=0.6, zorder=3, label="RMSD ≤ 2.0 Å success")
for xi, yi in zip(x, rmsd_ok):
    ax.text(xi, yi + 1.5, f"{yi:.1f}%", ha="center", va="bottom",
            fontsize=10.5, fontweight="bold", color=colors3[xi])

ax2 = ax.twinx()
nac_x = [xi for xi, v in zip(x, nac_ok) if v is not None]
nac_v = [v for v in nac_ok if v is not None]
ax2.plot(nac_x, nac_v, "-o", color="#185FA5", linewidth=2.2, markersize=7,
         zorder=4, label="NAC success rate")
for xi, vi in zip(nac_x, nac_v):
    ax2.annotate(f"{vi:.1f}%", xy=(xi, vi), xytext=(0, 10), textcoords="offset points",
                 ha="center", fontsize=10, color="#185FA5", fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel("RMSD ≤ 2.0 Å success rate (%)")
ax.set_ylim(0, 90)
ax2.set_ylabel("NAC success rate (%)", color="#185FA5")
ax2.set_ylim(0, 90)
ax2.tick_params(axis="y", colors="#185FA5")
ax2.spines[["top"]].set_visible(False)
ax.spines[["top"]].set_visible(False)
ax.grid(axis="y", color="#E5E7EB", linewidth=0.8, zorder=0)

h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, frameon=False, loc="upper left", fontsize=9.5)

fig.tight_layout()
fig.savefig(f"{OUT}/fig3_covalent.png", bbox_inches="tight")
plt.close(fig)

print("Figures written to", OUT)
import os
for f in sorted(os.listdir(OUT)):
    print(" ", f, os.path.getsize(os.path.join(OUT, f)), "bytes")
