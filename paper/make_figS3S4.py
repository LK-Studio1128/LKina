#!/usr/bin/env python3
"""Generate supplementary figures S3 (C3 deep-pocket ablation) and
S4 (metallocomplex full-pool) for the LKina manuscript.

S3: grouped bar of REACTIVE_DIST per preset for single / c3 / c3b, with the
    2.6 A convergence threshold dashed and per-variant convergence counts.
S4: per-system top-1 RMSD, LKina metal-as-ligand vs LKina AD4 standard,
    grouped by metal family; Vina parse-failure status annotated.

Output: figures/figS3_c3_ablation.{png,pdf}, figures/figS4_metallocomplex.{png,pdf}
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "legend.fontsize": 8,
    "figure.dpi": 600, "savefig.dpi": 600,
    "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
})

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
BENCH = os.path.join(HERE, "..", "byi", "LKina", "benchmarks")

VARIANT_COLORS = {"single": "#8b0000", "c3": "#b0b0b0", "c3b": "#2f5fa3"}
THRESH = 2.6  # A, deep-pocket convergence


def load_c3():
    d = json.load(open(os.path.join(BENCH, "c3_ablation_results.json")))
    presets = []
    dists = {"single": [], "c3": [], "c3b": []}
    conv = {"single": 0, "c3": 0, "c3b": 0}
    for preset, row in sorted(d["results"].items()):
        presets.append(preset)
        for v in ("single", "c3", "c3b"):
            e = row["variants"][v]
            dd = e.get("REACTIVE_DIST")
            dists[v].append(dd if dd is not None else 20.0)  # missing -> off-scale
            if e.get("converged"):
                conv[v] += 1
    return presets, dists, conv


def load_meta():
    rows = []
    for mc in ("PT", "PD", "RU", "OS", "RE"):
        f = os.path.join(BENCH, "metallocomplex_results", f"metallocomplex_{mc.lower()}.json")
        if not os.path.exists(f):
            continue
        d = json.load(open(f))
        for r in d.get("results", []):
            if r.get("status") != "ok":
                continue
            lk, ad = r.get("lkina_metal", {}), r.get("ad4_std", {})
            rows.append((mc, r["id"].split("_")[0], r.get("lig"),
                         lk.get("rmsd"), lk.get("metal_geom"), ad.get("rmsd")))
    return rows


def fig_s3(presets, dists, conv):
    x = np.arange(len(presets))
    w = 0.26
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    for i, v in enumerate(("single", "c3", "c3b")):
        bars = ax.bar(x + (i - 1) * w, dists[v], w, color=VARIANT_COLORS[v],
                      label=f"{v}  ({conv[v]}/{len(presets)} converged)",
                      edgecolor="#333333", linewidth=0.4)
        # annotate values on bars (skip off-scale)
        for xi, val in zip(x + (i - 1) * w, dists[v]):
            if val < 19:
                ax.text(xi, val + 0.15, f"{val:.1f}", ha="center", va="bottom", fontsize=6.5)
    ax.axhline(THRESH, color="#333333", ls="--", lw=0.9)
    ax.text(1, THRESH + 0.12, "convergence threshold 2.6 Å",
            ha="center", va="bottom", fontsize=7, style="italic")
    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("_", "\n") for p in presets], fontsize=7.5)
    ax.set_ylabel("Reactive distance d(SG–nucleophile) (Å)")
    ax.set_ylim(0, 15.5)
    ax.legend(frameon=False, ncol=1, loc="upper right", bbox_to_anchor=(0.99, 0.99),
              fontsize=7.2, handlelength=1.2, labelspacing=0.4)
    ax.set_title("C3 deep-pocket ablation — single-stage vs C3 vs C3b", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"figS3_c3_ablation.{ext}"))
    plt.close(fig)
    print("saved figS3_c3_ablation.{png,pdf}")


def fig_s4(rows):
    # order families, keep systems in pool order
    fam_order = ["PT", "PD", "RU", "OS", "RE"]
    fams = [f for f in fam_order if any(r[0] == f for r in rows)]
    labels, lk_rmsd, ad_rmsd = [], [], []
    ticks = []
    n = 0
    for fam in fams:
        sub = [r for r in rows if r[0] == fam]
        for r in sub:
            labels.append(f"{r[1]}\n{fam}·{r[2]}")
            lk_rmsd.append(r[3])
            ad_rmsd.append(r[5])
            ticks.append(n)
            n += 1
    x = np.arange(n)
    w = 0.38
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.bar(x - w / 2, lk_rmsd, w, color="#8b0000",
           label="LKina metal-as-ligand (20/20 docked)", edgecolor="#333333", lw=0.4)
    ax.bar(x + w / 2, ad_rmsd, w, color="#d9d2c5",
           label="LKina AD4 standard (--no_auto_metal)", edgecolor="#333333", lw=0.4)
    ax.axhline(2.0, color="#333333", ls="--", lw=0.9)
    ax.text(n - 0.4, 2.1, "redock success ≤ 2.0 Å", ha="right", va="bottom",
            fontsize=7, style="italic")
    for xi, v in zip(x, lk_rmsd):
        if v is not None:
            ax.text(xi - w / 2, v + 0.15, f"{v:.1f}", ha="center", va="bottom", fontsize=6)
    for xi, v in zip(x, ad_rmsd):
        if v is not None:
            ax.text(xi + w / 2, v + 0.15, f"{v:.1f}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Top-1 RMSD vs crystal pose (Å)")
    ax.set_ylim(0, 12.5)
    ax.legend(frameon=False, ncol=1, loc="upper left", fontsize=7.5)
    ax.set_title("Metal-as-ligand redocking — full Pt/Pd/Ru/Os/Re pool (n=20); "
                 "Vina 1.2.7 fails 20/20 at parse time", fontsize=9.5)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"figS4_metallocomplex.{ext}"))
    plt.close(fig)
    print("saved figS4_metallocomplex.{png,pdf}")


if __name__ == "__main__":
    presets, dists, conv = load_c3()
    fig_s3(presets, dists, conv)
    rows = load_meta()
    print(f"S4 data rows: {len(rows)}")
    fig_s4(rows)
