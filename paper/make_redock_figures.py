#!/usr/bin/env python3
"""Regenerate Figure 1 / Figure 2 (redocking benchmark) from archived measured data.

Fig 1: per-metal top-1 RMSD success + donor-metal <=3A rates (3 engines)
Fig 2: per-system best-pose donor-metal distance, 104 systems, 3 engines

All numbers recomputed from benchmarks/redock_benchmark/results/*.json and
docking/*.pdbqt so figures match the paper tables exactly.
"""
import json, glob, os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks/redock_benchmark"
RES  = os.path.join(BASE, "results")
DOCK = os.path.join(BASE, "docking")
OUT  = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/figures"

ENGINES = ["lkina_metal", "ad4_std", "vina127"]
LABELS  = {"lkina_metal": "LKina AD4 + metal mode",
           "ad4_std":     "LKina AD4 (--no_auto_metal)",
           "vina127":     "Vina 1.2.7"}
COLORS  = {"lkina_metal": "#C0392B",   # LKina red
           "ad4_std":     "#5A8F5A",   # green (C tier, darkened)
           "vina127":     "#6B7FA3"}   # grey-blue (Vina)
DONOR_TYPES = {"N", "NA", "NS", "OA", "OS", "S", "SA"}  # ligand-side donors incl. AD4 types

# ---------------------------------------------------------------- load data
files = {"ZN": f"{RES}/redock_zn_120.json",
         "FE": f"{RES}/redock_fe_120.json",
         "CU": f"{RES}/redock_cu_120.json"}
data, meta_all = {}, {}
for metal, path in files.items():
    data[metal] = json.load(open(path))
    for r in data[metal]["results"]:
        mp = f"{RES}/{r['id']}_meta.json"
        if os.path.exists(mp):
            meta_all[r["id"]] = json.load(open(mp))

def parse_pose_atoms(path):
    """Return list of (x,y,z,ad4type) from first MODEL of a docked pdbqt.

    Column-exact [77:79] parsing breaks when the ligand name is 4 characters
    (e.g. 0ZB) — fall back to whitespace tokenization for the trailing type.
    """
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL") or line.rstrip() == "END":
                break
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                t = line[77:79].strip()
                if not t:  # column-shifted layout → take the last token
                    toks = line.split()
                    t = toks[-1] if toks else ""
                atoms.append((x, y, z, t))
    return atoms

def pose_donor_metal_min(tag, engine):
    """min over pose N/O/S heavy atoms of distance to nearest receptor metal."""
    metals = []
    for mid, m in meta_all.items():
        if m["pdb"] == tag.split("_")[0] and f"{m['metal']}" == tag.split("_")[1]:
            metals = [tuple(c) for c in m["metals"]]
            break
    if not metals:
        return None
    suffix = {"lkina_metal": None, "ad4_std": "ad4", "vina127": "vina"}[engine]
    if suffix is None:  # lkina metal-mode output uses the mode token (zn/fe3/cu2_jt)
        mode = {"ZN": "zn", "FE": "fe3", "CU": "cu2_jt"}[tag.split("_")[1]]
        dpath = f"{DOCK}/{tag}_lkina_{mode}.pdbqt"
    else:
        dpath = f"{DOCK}/{tag}_{suffix}.pdbqt"
    if not os.path.exists(dpath):
        return None
    try:
        atoms = parse_pose_atoms(dpath)
    except Exception:
        return None
    donors = [(x, y, z) for (x, y, z, t) in atoms if t in DONOR_TYPES]
    if not donors:
        return None
    return min(math.dist(d, mt) for d in donors for mt in metals)

# ---------------------------------------------------------------- recompute
per_system = {e: [] for e in ENGINES}          # list of (id, dm_dist)
for metal in files:
    for r in data[metal]["results"]:
        tag = r["id"]
        # only systems where all three engines docked a pose with finite rmsd/metadata exist
        ok = all(r[e].get("rmsd") is not None for e in ENGINES)
        if not ok or tag not in meta_all:
            continue
        for e in ENGINES:
            d = pose_donor_metal_min(tag, e)
            if d is not None:
                per_system[e].append((tag, d))

summary_agg = {e: {
    "n": len(v),
    "mean": float(np.mean([d for _, d in v])),
    "median": float(np.median([d for _, d in v])),
    "le3": int(sum(1 for _, d in v if d <= 3.0)),
} for e, v in per_system.items()}
print(json.dumps(summary_agg, indent=1))

# cross-check against archived summary
arch = json.load(open(f"{RES}/donor_metal_distance_summary.json"))
for e in ENGINES:
    a, c = arch[e], summary_agg[e]
    print(f"{e}: recomputed n={c['n']} median={c['median']:.2f} le3={c['le3']} | "
          f"archived n={a['n']} median={a['median']} le3={a['le_3A']}")

# ---------------------------------------------------------------- styling
plt.rcParams.update({
    "font.size": 10, "axes.linewidth": 0.9,
    "savefig.dpi": 600,
})
plt.rcParams["font.family"] = ["Helvetica", "Arial", "DejaVu Sans"]

# ================================================================ FIGURE 1
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))

metals  = ["ZN", "FE", "CU"]
mlabels = ["Zn\u00b2\u207a\n(n=22)", "Fe\u00b3\u207a\n(n=47)", "Cu\u00b2\u207a\n(n=35)"]
ax = axes[0]
x = np.arange(3); w = 0.26
rmsd_tbl = json.load(open(f"{RES}/redock_summary.json"))
for i, e in enumerate(ENGINES):
    vals = [rmsd_tbl[m][e]["rate_pct"] for m in metals]
    xs = x + (i - 1) * w
    bars = ax.bar(xs, vals, w * 0.92, color=COLORS[e], label=LABELS[e],
                  edgecolor="white", linewidth=0.5, zorder=3)
    for xi, yi in zip(xs, vals):
        ax.text(xi, yi + 0.8, f"{yi:.0f}", ha="center", va="bottom",
                fontsize=8.5, fontweight="bold", color=COLORS[e])
ax.set_xticks(x); ax.set_xticklabels(mlabels, fontsize=9.5)
ax.set_ylabel("Top-1 RMSD \u2264 2.0 \u00c5 (%)")
ax.set_ylim(0, 42)
ax.grid(axis="y", color="#E5E7EB", linewidth=0.8, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="upper right")
ax.set_title("(A) Global-pose success", fontsize=10.5)

ax = axes[1]
for i, e in enumerate(ENGINES):
    vals = [100.0 * sum(1 for _, d in per_system[e] if d <= 3.0 and _.split('_')[1] == m) /
            max(1, sum(1 for t, _ in per_system[e] if t.split('_')[1] == m)) for m in metals]
    xs = x + (i - 1) * w
    bars = ax.bar(xs, vals, w * 0.92, color=COLORS[e], label=LABELS[e],
                  edgecolor="white", linewidth=0.5, zorder=3)
    for xi, yi in zip(xs, vals):
        ax.text(xi, yi + 1.2, f"{yi:.0f}", ha="center", va="bottom",
                fontsize=8.5, fontweight="bold", color=COLORS[e])
ax.set_xticks(x); ax.set_xticklabels(mlabels, fontsize=9.5)
ax.set_ylabel("Pose donor within 3.0 \u00c5 of metal (%)")
ax.set_ylim(0, 108)
ax.grid(axis="y", color="#E5E7EB", linewidth=0.8, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("(B) Coordination-competent poses", fontsize=10.5)

fig.tight_layout()
fig.savefig(f"{OUT}/fig6_redock_benchmark.png", bbox_inches="tight")
fig.savefig(f"{OUT}/fig6_redock_benchmark.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/fig6_redock_benchmark.tif", bbox_inches="tight")
plt.close(fig)

# ================================================================ FIGURE 2
fig, ax = plt.subplots(figsize=(8.6, 4.2))
order = sorted({t for t, _ in per_system["lkina_metal"]},
               key=lambda s: per_system["lkina_metal"][[t for t,_ in per_system["lkina_metal"]].index(s)][1])
xs = np.arange(len(order))
for i, e in enumerate(ENGINES):
    lut = dict(per_system[e])
    ys = [lut.get(t, np.nan) for t in order]
    ms = 16 if e != "lkina_metal" else 20
    ax.scatter(xs, ys, s=ms, color=COLORS[e], alpha=0.75 if e != "lkina_metal" else 0.95,
               linewidths=0, zorder=3, label=LABELS[e])

meds = {e: float(np.median([d for _, d in v])) for e, v in per_system.items()}
for j, e in enumerate(ENGINES):
    ax.axhline(meds[e], color=COLORS[e], linestyle="--", linewidth=1.1, alpha=0.8, zorder=2)
    ax.text(len(order) + 2, meds[e], f"median {meds[e]:.2f} \u00c5",
            color=COLORS[e], fontsize=8.5, va="center", fontweight="bold")

ax.axhspan(0, 3.0, color="#FCEBEB", zorder=0)
ax.text(-12, 1.5, "coordination sphere \u2264 3.0 \u00c5", fontsize=8.5,
        color="#C0392B", ha="left", va="center", fontweight="bold")
ax.set_xlim(-4, len(order) + 1)
ax.set_xticks([0, len(order) // 2, len(order) - 1])
ax.set_xticklabels([order[0].split("_")[0], order[len(order) // 2].split("_")[0], order[-1].split("_")[0]])
ax.set_xlabel("Metalloprotein system (104 total: Zn 22 / Fe 47 / Cu 35)")
ax.set_ylabel("Best-pose donor\u2013metal min distance (\u00c5)")
ax.set_ylim(0, 14)
ax.grid(axis="y", color="#EEF0F3", linewidth=0.6, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=8.5, loc="upper left", markerscale=2.2)
fig.tight_layout()
fig.savefig(f"{OUT}/fig7_donor_metal_distance.png", bbox_inches="tight")
fig.savefig(f"{OUT}/fig7_donor_metal_distance.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/fig7_donor_metal_distance.tif", bbox_inches="tight")
plt.close(fig)
print("figures written:", OUT)
