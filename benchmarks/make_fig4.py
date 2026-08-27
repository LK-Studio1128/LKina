#!/usr/bin/env python3
"""Figure 4: LKina vs Vina 1.2.7 metal-type coverage & coordination recovery."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import json

OUT = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/figures"
data = json.load(open("/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks/metal_coverage_results_v3.json"))

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.edgecolor": "#444", "axes.linewidth": 0.9,
    "legend.fontsize": 10, "xtick.labelsize": 9.5, "ytick.labelsize": 10,
    "figure.dpi": 200,
})

lk_ok = [r for r in data if r["lkina_rc"] == 0]
vina_ok = [r for r in data if r["vina_rc"] == 0]
n = len(data)

# ---- Panel A: coverage bar ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2), gridspec_kw={"width_ratios": [1, 1.6]})

bars = ax1.bar(["LKina", "Vina 1.2.7"],
               [len(lk_ok), len(vina_ok)],
               color=["#C0392B", "#9AA5B1"], width=0.55, edgecolor="white")
ax1.text(0, len(lk_ok) + 1, f"{len(lk_ok)}/63", ha="center", fontweight="bold", color="#C0392B")
ax1.text(1, len(vina_ok) + 1, f"{len(vina_ok)}/63", ha="center", fontweight="bold", color="#666")
ax1.set_ylim(0, 75)
ax1.set_ylabel("Metal modes dockable")
ax1.set_title("Metal-mode coverage (63 modes)")
ax1.spines[["top", "right"]].set_visible(False)
ax1.grid(axis="y", color="#E5E7EB", linewidth=0.8, zorder=0)
ax1.annotate("47 modes fail:\n\"Atom type Pt is not\na valid AutoDock type\"",
             xy=(1, len(vina_ok)), xytext=(1.32, 40), fontsize=8.5, color="#791F1F",
             ha="left", va="center", fontweight="bold",
             arrowprops=dict(arrowstyle="->", color="#791F1F", lw=0.9))
ax1.set_xlim(-0.6, 2.4)

# ---- Panel B: coordination distance error ----
modes, lk_err, vina_err = [], [], []
for r in sorted(data, key=lambda x: x["mode"]):
    if r["lkina_rc"] != 0: continue
    modes.append(r["mode"])
    lk_err.append(abs((r["d_lig"] or 99) - r["d0"]))
    vina_err.append(abs((r["vina_d"] or 99) - r["d0"]) if r["vina_rc"] == 0 else None)

x = np.arange(len(modes))
ax2.bar(x - 0.2, lk_err, 0.4, color="#C0392B", label="LKina (AD4+metal)", edgecolor="white")
vina_vals = [v if v is not None else 0 for v in vina_err]
vina_mask = [v is not None for v in vina_err]
ax2.bar(x + 0.2, [vina_vals[i] if vina_mask[i] else 0 for i in range(len(modes))],
        0.4, color="#9AA5B1", label="Vina 1.2.7 (dockable only)", edgecolor="white")
# mark vina-failed modes with hatch overlay
ax2.bar(x + 0.2, [0.3 if not vina_mask[i] else 0 for i in range(len(modes))],
        0.4, color="none", edgecolor="#666", hatch="//", linewidth=0.6)

ax2.set_xticks(x)
ax2.set_xticklabels(modes, rotation=90, fontsize=7)
ax2.set_ylabel("|d(M–donor) − d₀| (Å)")
ax2.set_title("Coordination-distance error (lower = better)")
ax2.axhline(0.5, color="#C0392B", linestyle=":", linewidth=1.0)
ax2.text(len(modes)*0.98, 0.55, "0.5 Å", color="#C0392B", fontsize=8, ha="right")
ax2.set_ylim(0, 8)
ax2.spines[["top", "right"]].set_visible(False)
ax2.grid(axis="y", color="#E5E7EB", linewidth=0.8, zorder=0)
ax2.legend(frameon=False, loc="upper right", fontsize=9)
# hatched = vina cannot dock
ax2.text(0.02, 0.97, "hatched: Vina cannot parse metal type",
         transform=ax2.transAxes, fontsize=8, color="#666", va="top")

fig.tight_layout()
fig.savefig(f"{OUT}/fig4_metal_coverage.png", bbox_inches="tight")
plt.close(fig)
print("fig4 saved")

# ---- summary stats ----
errs = lk_err
import statistics
print(f"LKina: {sum(1 for e in errs if e < 0.5)}/63 err<0.5A, {sum(1 for e in errs if e < 1.0)}/63 err<1.0A")
print(f"mean err {statistics.mean(errs):.2f} median {statistics.median(errs):.2f}")
# only modes where both ran
both = [(a, b) for a, b in zip(lk_err, vina_err) if b is not None]
lk_b = [a for a, b in both]; vina_b = [b for a, b in both]
print(f"both-ran {len(both)} modes: LKina mean err {statistics.mean(lk_b):.2f}, Vina mean err {statistics.mean(vina_b):.2f}")
print(f"LKina better: {sum(1 for a,b in both if a < b)}/{len(both)}")
