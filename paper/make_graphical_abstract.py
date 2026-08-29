#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Graphical Abstract (TOC) for the LKina paper.

Concept layout (single visual, JCIM style):
  TOP    : AutoDock Vina 1.2.7 -> LKina engine (inheritance arrow)
  LEFT   : metal-coordination docking  (Zn2+ + pseudoatoms, 4JC Zn-NA 2.13 A)
  RIGHT  : covalent-reactive docking   (Cys-SG -> electrophile, NAC angle)
  BOTTOM : inline AutoGrid (no external autogrid4) + backward compatible
  STRIP  : Vina 1.2.7 34/110  ->  LKina 110/110
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Ellipse, Arc
import numpy as np, os

FIG = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/figures"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "figure.dpi": 600, "savefig.dpi": 600,
    "axes.linewidth": 0.8,
})

RED   = "#8b0000"
RED2  = "#c0392b"
CREAM = "#efe9dc"
GOLD  = "#c9a227"
GREY  = "#555555"
LGRY  = "#cfcfcf"
ZNC   = "#3f6fb5"
BLK   = "#1a1a1a"

fig, ax = plt.subplots(figsize=(6.5, 3.3))
ax.set_xlim(0, 13); ax.set_ylim(0, 6.6); ax.axis("off")

# ---------------- top: Vina -> LKina ----------------
ax.add_patch(FancyBboxPatch((4.7, 5.55), 3.6, 0.78, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc="#f7f5f0", ec=LGRY, lw=0.8, zorder=2))
ax.text(6.5, 5.94, "AutoDock Vina 1.2.7  (32-type XS)", ha="center", va="center",
        fontsize=8.5, color=GREY)
ax.add_patch(FancyArrowPatch((6.5, 5.45), (6.5, 4.92), arrowstyle="-|>", mutation_scale=16,
                             color=RED, lw=2.0, zorder=3))
ax.text(7.0, 5.18, "extends", fontsize=7, color=RED, style="italic")

ax.add_patch(FancyBboxPatch((4.1, 4.02), 4.8, 0.82, boxstyle="round,pad=0.05,rounding_size=0.14",
                            fc=RED, ec=RED, lw=1.2, zorder=2))
ax.text(6.5, 4.43, "LKina", ha="center", va="center", fontsize=13, fontweight="bold", color="white")
ax.text(6.5, 4.14, "113 types  ·  metal-aware  ·  covalent-reactive", ha="center", va="center",
        fontsize=7.2, color="#f3d9d9")

# ---------------- LEFT: metal coordination ----------------
ax.add_patch(FancyBboxPatch((0.3, 1.45), 3.7, 2.4, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc=CREAM, ec=RED2, lw=1.0, zorder=1))
ax.text(2.15, 3.62, "Metal-coordination docking", ha="center", fontsize=8.5,
        fontweight="bold", color=RED)
# pocket + Zn + ligand (raised, slightly smaller, clear of title and text)
pocket = Ellipse((1.15, 3.04), 1.0, 0.82, angle=25, fc="#e7e2d6", ec=GREY, lw=0.7, zorder=2)
ax.add_patch(pocket)
ax.add_patch(Circle((1.15, 3.04), 0.18, fc=ZNC, ec="white", lw=1.0, zorder=4))
ax.text(1.15, 3.04, "Zn", ha="center", va="center", color="white", fontsize=6.5, fontweight="bold", zorder=5)
# pseudoatoms (small dashed circles at coordination sites)
for ang in (20, 140, 260):
    a = np.deg2rad(ang)
    ax.add_patch(Circle((1.15 + 0.43*np.cos(a), 3.04 + 0.43*np.sin(a)), 0.068,
                        fc="none", ec=RED2, lw=1.0, ls=(0,(2,1)), zorder=4))
# ligand N donor approaching (Zn-NA 2.13 A)
ax.add_patch(Circle((1.77, 3.27), 0.12, fc=GOLD, ec="white", lw=0.8, zorder=4))
ax.text(1.77, 3.27, "N", ha="center", va="center", color="white", fontsize=6, zorder=5)
ax.plot([1.32, 1.66], [3.08, 3.24], color=ZNC, lw=1.4, ls="--", zorder=3)
ax.text(1.52, 2.78, "2.13 Å", fontsize=6.5, color=ZNC, rotation=18)
# text stack (kept inside box, clear margin from bottom edge)
ax.text(2.15, 2.38, "TZ / SQ / MH / JT pseudoatoms", ha="center", fontsize=7.2, color=BLK)
ax.text(2.15, 2.15, "BVS oxidation-state inference", ha="center", fontsize=7.2, color=BLK)
ax.text(2.15, 1.92, "4JC Zn–NA at 2.13 Å", ha="center", fontsize=7.0, color=GREY)
ax.text(2.15, 1.70, "(Vina 2.23 Å, no coord)", ha="center", fontsize=6.3, color=GREY)

# ---------------- RIGHT: covalent reactive ----------------
ax.add_patch(FancyBboxPatch((9.0, 1.45), 3.7, 2.4, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc=CREAM, ec=RED2, lw=1.0, zorder=1))
ax.text(10.85, 3.62, "Covalent-reactive docking", ha="center", fontsize=8.5,
        fontweight="bold", color=RED)
# Cys-SG -> electrophile with NAC angle
pocket2 = Ellipse((11.5, 3.00), 1.1, 0.86, angle=-20, fc="#e7e2d6", ec=GREY, lw=0.7, zorder=2)
ax.add_patch(pocket2)
ax.add_patch(Circle((11.15, 3.16), 0.14, fc="#e0a106", ec="white", lw=0.8, zorder=4))
ax.text(11.15, 3.16, "SG", ha="center", va="center", color="white", fontsize=5.5, zorder=5)
ax.add_patch(Circle((11.74, 2.86), 0.12, fc=RED2, ec="white", lw=0.8, zorder=4))
ax.text(11.74, 2.86, "E$^{+}$", ha="center", va="center", color="white", fontsize=6, zorder=5)
ax.plot([11.28, 11.63], [3.11, 2.91], color=RED2, lw=1.6, zorder=3)
# NAC angle arc
ax.add_patch(Arc((11.15, 3.16), 0.70, 0.70, angle=0, theta1=-42, theta2=0,
                 color=GOLD, lw=1.2, zorder=3))
ax.text(11.60, 3.34, "NAC θ", fontsize=6.3, color=GOLD)
# text stack
ax.text(10.85, 2.38, "NAC detection (d < 3.0 Å, angle)", ha="center", fontsize=7.2, color=BLK)
ax.text(10.85, 2.15, "P1–P4 · C3 two-step · 6 presets", ha="center", fontsize=7.2, color=BLK)
ax.text(10.85, 1.92, "end-to-end, NAC discriminates", ha="center", fontsize=6.6, color=GREY)
ax.text(10.85, 1.71, "covalent vs non-covalent", ha="center", fontsize=6.3, color=GREY)

# ---------------- connectors from LKina to two wings ----------------
ax.add_patch(FancyArrowPatch((4.9, 4.10), (3.95, 3.42), arrowstyle="-|>", mutation_scale=13,
                             color=RED, lw=1.4, zorder=3))
ax.add_patch(FancyArrowPatch((8.1, 4.10), (9.05, 3.42), arrowstyle="-|>", mutation_scale=13,
                             color=RED, lw=1.4, zorder=3))

# ---------------- BOTTOM bar: inline grid + backward compat ----------------
ax.add_patch(FancyBboxPatch((0.3, 0.55), 6.2, 0.8, boxstyle="round,pad=0.04,rounding_size=0.1",
                            fc="#f4f2ee", ec=LGRY, lw=0.8, zorder=1))
ax.text(3.4, 1.12, "Inline AutoGrid 4.2  —  no external autogrid4", ha="center",
        fontsize=7.6, color=BLK)
ax.text(3.4, 0.82, "108 AG4-format maps incl. TZ/SQ/MH/JT", ha="center", fontsize=6.8, color=GREY)
ax.add_patch(FancyBboxPatch((6.8, 0.55), 5.9, 0.8, boxstyle="round,pad=0.04,rounding_size=0.1",
                            fc="#f4f2ee", ec=LGRY, lw=0.8, zorder=1))
ax.text(9.75, 1.12, "Backward compatible — Vina 1.2.7 energy identical", ha="center",
        fontsize=7.6, color=BLK)
ax.text(9.75, 0.82, "1HVR −14.54 · 3PTB −6.202 kcal/mol (mode-1)", ha="center", fontsize=6.8, color=GREY)

# ---------------- STRIP: Vina vs LKina coverage ----------------
ax.add_patch(FancyBboxPatch((0.3, 0.02), 12.4, 0.4, boxstyle="round,pad=0.03,rounding_size=0.06",
                            fc=RED, ec=RED, zorder=1))
ax.text(6.5, 0.22, "Metal-mode coverage:  Vina 1.2.7  34/110  →  LKina  110/110",
        ha="center", va="center", fontsize=8.2, color="white", fontweight="bold")

fig.tight_layout(pad=0.2)
for ext in ("png", "pdf", "tif"):
    fig.savefig(os.path.join(FIG, f"graphical_abstract.{ext}"), bbox_inches="tight")
plt.close(fig)
print("saved graphical_abstract.png/.pdf")
