#!/usr/bin/env python3
"""Generate all supplementary-test figures (fig S5-S10 series).

S5  scoring function comparison (energy + runtime + Zn distance)
S6  exhaustiveness convergence (energy + success vs exh, runtime scaling)
S7  seed variance (box/violin of best energy per seed)
S8  metal_soft_weight sweep (Zn-ligand distance curve)
S9  metal_bias strength sweep (Zn-ligand distance + energy)
S10 hotspot cases (6OIM covalent modes; 3HS4 zinc modes) summary

Inputs: scoring_comparison_results.json, parameter_sweep_results.json,
        hotspot_results.json
Output: figures/figS5..S10 (.png + .pdf + .tif), 600 dpi
"""
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BENCH = Path(__file__).parent
FIGS = BENCH.parent / "figures"
FIGS.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
    "figure.dpi": 600, "savefig.dpi": 600, "axes.spines.top": False,
    "axes.spines.right": False,
})
C = {"vina": "#4C72B0", "vinardo": "#DD8452", "ad4": "#55A868",
     "blue": "#4C72B0", "orange": "#DD8452", "green": "#55A868",
     "red": "#C44E52", "purple": "#8172B3", "gray": "#8C8C8C"}


def save(fig, name):
    fig.savefig(FIGS / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIGS / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS / f"{name}.tif", bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


# ---------------------------------------------------------------- figS5 scoring
def figS5():
    data = json.load(open(BENCH / "scoring_comparison_results.json"))
    ids = ["4JC", "3HS4", "6OIM"]
    scorings = ["vina", "vinardo", "ad4"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    x = np.arange(len(ids))
    w = 0.25

    best = {s: [] for s in scorings}
    rt = {s: [] for s in scorings}
    zd = {s: [] for s in scorings}
    lookup = {(d["id"], d["scoring"]): d for d in data}
    for s in scorings:
        for i in ids:
            d = lookup.get((i, s), {})
            best[s].append(d.get("best") or np.nan)
            rt[s].append(d.get("runtime_s") or np.nan)
            v = d.get("min_metal_dist_A")
            zd[s].append(v if v is not None else np.nan)

    ax = axes[0]
    for k, s in enumerate(scorings):
        ax.bar(x + (k-1)*w, best[s], w, label=s, color=C[s])
    ax.set_xticks(x); ax.set_xticklabels(ids)
    ax.set_ylabel("Best affinity (kcal/mol)")
    ax.set_title("(a) Best score")

    ax = axes[1]
    for k, s in enumerate(scorings):
        ax.bar(x + (k-1)*w, rt[s], w, label=s, color=C[s])
    ax.set_xticks(x); ax.set_xticklabels(ids)
    ax.set_ylabel("Runtime (s, exh=8, 1 mode)")
    ax.set_title("(b) Runtime")

    ax = axes[2]
    for k, s in enumerate(scorings):
        ax.bar(x + (k-1)*w, zd[s], w, label=s, color=C[s])
    ax.set_xticks(x); ax.set_xticklabels(ids)
    ax.set_ylabel("Min ligand-to-Zn distance (Å)")
    ax.set_title("(c) Metal proximity")
    ax.axhline(2.1, ls="--", lw=0.8, color=C["red"])
    ax.text(2.1, 2.2, "Zn–N ref", fontsize=8, color=C["red"])

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout()
    save(fig, "figS5_scoring_comparison")


# ---------------------------------------------------------------- figS6 exh
def figS6():
    data = json.load(open(BENCH / "parameter_sweep_results.json"))["A_exhaustiveness"]
    exhs = sorted({d["exhaustiveness"] for d in data})
    meanE, allE, meanT = [], [], []
    for e in exhs:
        rows = [d for d in data if d["exhaustiveness"] == e and d.get("best") is not None]
        allE.append([d["best"] for d in rows])
        meanE.append(np.mean([d["best"] for d in rows]) if rows else np.nan)
        meanT.append(np.mean([d["runtime_s"] for d in rows if d.get("runtime_s")]) if rows else np.nan)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    ax = axes[0]
    bp = ax.boxplot(allE, positions=range(len(exhs)), widths=0.5, patch_artist=True,
                    medianprops=dict(color="k"))
    for p in bp["boxes"]:
        p.set_facecolor(C["blue"]); p.set_alpha(0.6)
    ax.set_xticklabels(exhs)
    ax.set_xlabel("exhaustiveness"); ax.set_ylabel("Best affinity (kcal/mol)")
    ax.set_title("(a) Score convergence (3HS4 Zn, 3 seeds)")

    ax = axes[1]
    sd = [np.std(v) if len(v) > 1 else 0 for v in allE]
    ax.plot(exhs, sd, "o-", color=C["red"])
    ax.set_xlabel("exhaustiveness"); ax.set_ylabel("σ of best affinity (kcal/mol)")
    ax.set_title("(b) Seed-to-seed spread")

    ax = axes[2]
    ax.plot(exhs, meanT, "s-", color=C["green"])
    ax.set_xlabel("exhaustiveness"); ax.set_ylabel("Runtime (s)")
    ax.set_title("(c) Runtime scaling")
    fig.tight_layout()
    save(fig, "figS6_exhaustiveness")


# ---------------------------------------------------------------- figS7 seeds
def figS7():
    data = json.load(open(BENCH / "parameter_sweep_results.json"))["B_seed_variance"]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    systems = ["3HS4", "6OIM"]
    for k, sys in enumerate(systems):
        rows = [(d["seed"], d["best"]) for d in data
                if d.get("system") == sys and d.get("best") is not None]
        rows.sort()
        seeds = [r[0] for r in rows]; vals = [r[1] for r in rows]
        ax.plot(seeds, vals, "o-", color=list(C.values())[k], label=sys)
        for s, v in zip(seeds, vals):
            ax.annotate(f"{v:.2f}", (s, v), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax.set_xlabel("random seed"); ax.set_ylabel("Best affinity (kcal/mol)")
    ax.set_title("Seed variance at exhaustiveness 8")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "figS7_seed_variance")


# ---------------------------------------------------------------- figS8 soft weight
def figS8():
    data = json.load(open(BENCH / "parameter_sweep_results.json"))["C_soft_weight"]
    ws, dist, E = [], [], []
    for d in sorted(data, key=lambda x: x["soft_weight"]):
        ws.append(d["soft_weight"])
        dist.append(d.get("min_zn_dist") or np.nan)
        E.append(d.get("best") or np.nan)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(ws, dist, "o-", color=C["purple"], label="min ligand-to-Zn distance")
    ax.set_xlabel("--metal_soft_weight"); ax.set_ylabel("Distance (Å)", color=C["purple"])
    ax2 = ax.twinx()
    ax2.plot(ws, E, "s--", color=C["orange"], label="best affinity")
    ax2.set_ylabel("Affinity (kcal/mol)", color=C["orange"])
    ax2.spines["right"].set_visible(True)
    ax.set_title("Metal soft-constraint gradient weight sweep (3HS4, zn)")
    fig.tight_layout()
    save(fig, "figS8_soft_weight")


# ---------------------------------------------------------------- figS9 bias
def figS9():
    data = json.load(open(BENCH / "parameter_sweep_results.json"))["D_metal_bias"]
    ss, dist, E = [], [], []
    for d in sorted(data, key=lambda x: x["bias_strength"]):
        ss.append(d["bias_strength"])
        dist.append(d.get("min_zn_dist") or np.nan)
        E.append(d.get("best") or np.nan)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(ss, dist, "o-", color=C["purple"], label="min ligand-to-Zn distance")
    ax.set_xlabel("--metal_bias_strength (kcal/mol well depth)")
    ax.set_ylabel("Distance (Å)", color=C["purple"])
    ax2 = ax.twinx()
    ax2.plot(ss, E, "s--", color=C["orange"], label="best affinity")
    ax2.set_ylabel("Affinity (kcal/mol)", color=C["orange"])
    ax2.spines["right"].set_visible(True)
    ax.set_title("Metal-bias (O5) attractor strength sweep (3HS4, zn)")
    fig.tight_layout()
    save(fig, "figS9_metal_bias")


# ---------------------------------------------------------------- figS10 hotspot
def figS10():
    hot = json.load(open(BENCH / "hotspot_results.json"))
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))

    # 6OIM: bars for best energy + line for warhead-SG distance
    rows = hot["6OIM_covalent"]
    labels = [r["tag"].replace("6OIM_", "") for r in rows]
    best = [r.get("best") or np.nan for r in rows]
    dsg = [r.get("warhead_C17_to_SG") or np.nan for r in rows]
    x = np.arange(len(rows))
    ax = axes[0]
    ax.bar(x - 0.18, best, 0.36, color=C["blue"], label="best affinity")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=12)
    ax.set_ylabel("Affinity (kcal/mol)", color=C["blue"])
    ax2 = ax.twinx()
    ax2.plot(x, dsg, "o--", color=C["red"], label="C17(warhead)–SG(C12) dist")
    ax2.set_ylabel("Warhead–Cys12 S distance (Å)", color=C["red"])
    ax2.spines["right"].set_visible(True)
    ax.set_title("(a) KRAS G12C (6OIM): covalent modes")
    ax.axhline(np.nanmin(best) if not all(math.isnan(b) for b in best) else 0, lw=0)

    rows = hot["3HS4_zinc"]
    labels = [r["tag"].replace("3HS4_", "") for r in rows]
    best = [r.get("best") or np.nan for r in rows]
    dz = [r.get("min_dist_to_Zn") or np.nan for r in rows]
    ds = [r.get("S_to_Zn_dist") or np.nan for r in rows]
    x = np.arange(len(rows))
    ax = axes[1]
    ax.bar(x - 0.18, best, 0.36, color=C["green"], label="best affinity")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=12)
    ax.set_ylabel("Affinity (kcal/mol)", color=C["green"])
    ax2 = ax.twinx()
    ax2.plot(x, dz, "o--", color=C["purple"], label="ligand–Zn min dist")
    ax2.plot(x, ds, "^:", color=C["red"], label="AZM S–Zn dist")
    ax2.set_ylabel("Distance to Zn (Å)")
    ax2.spines["right"].set_visible(True)
    ax2.legend(frameon=False, fontsize=8)
    ax.set_title("(b) CA II (3HS4): zinc-aware modes")

    fig.tight_layout()
    save(fig, "figS10_hotspot_cases")


if __name__ == "__main__":
    import sys
    only = sys.argv[1] if len(sys.argv) > 1 else None
    fns = {"S5": figS5, "S6": figS6, "S7": figS7, "S8": figS8, "S9": figS9, "S10": figS10}
    if only:
        fns[only]()
    else:
        for f in fns.values():
            try:
                f()
            except FileNotFoundError as e:
                print("skip:", e)
