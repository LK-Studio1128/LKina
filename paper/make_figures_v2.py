#!/usr/bin/env python3
"""Generate the supplementary measured-data figures for the LKina manuscript:
Fig 5  all-metal coverage (110 modes, LKina vs Vina 1.2.7)
Fig 6  coordination accuracy + pseudoatom geometry
Fig 7  BVS oxidation-state inference
Fig 8  covalent framework (P1-P4, C3, NAC, well scan)
Fig 9  inline map generation + metal-as-ligand geometry QC
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import json, os

BENCH = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks"
OUT   = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.edgecolor": "#444",
    "axes.linewidth": 0.9,
    "axes.titlesize": 12,
    "axes.labelsize": 11.5,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 200,
})
C_LKINA = "#C0392B"
C_VINA  = "#9AA5B1"
C_GREY  = "#B8BEC6"
C_GOLD  = "#E8A33D"
C_BLUE  = "#2E6FAD"

metal = json.load(open(os.path.join(BENCH, "metal_coverage_results_all.json")))
ff    = json.load(open(os.path.join(BENCH, "feature_family_results.json")))
cov   = json.load(open(os.path.join(BENCH, "covalent_full_results.json")))

# ---------------------------------------------------------------------------
# Figure 5 — all-metal coverage: LKina 110/110 vs Vina 1.2.7 34/110
# ---------------------------------------------------------------------------
n_ok = sum(1 for r in metal if r["lkina_rc"] == 0)
n_vok = sum(1 for r in metal if r["vina_rc"] == 0)
n_fail = sum(1 for r in metal if r["vina_rc"] != 0)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.4), gridspec_kw={"width_ratios": [1.25, 1]})

x = np.arange(2)
w = 0.55
ax1.bar(x, [n_ok, n_vok], w, color=[C_LKINA, C_VINA], edgecolor="#333", linewidth=0.6)
ax1.text(0, n_ok + 2, f"{n_ok}/110", ha="center", fontweight="bold", color=C_LKINA, fontsize=13)
ax1.text(1, n_vok + 2, f"{n_vok}/110", ha="center", fontweight="bold", color="#666", fontsize=13)
ax1.annotate("76 metal types\nrejected at PDBQT parse\n('Atom type X is not\na valid AutoDock type')",
             xy=(1, n_vok), xytext=(0.42, 52), fontsize=9, color="#333",
             arrowprops=dict(arrowstyle="->", color="#666", lw=1.0))
ax1.set_xticks(x); ax1.set_xticklabels(["LKina", "Vina 1.2.7"])
ax1.set_ylabel("Metal coordination modes docked successfully")
ax1.set_ylim(0, 125)
ax1.set_title("(a) Coverage: 110 metal modes", fontsize=11.5)
ax1.spines[["top", "right"]].set_visible(False)
ax1.grid(axis="y", color="#E5E7EB", linewidth=0.8, zorder=0)

# right panel: Vina-ok vs Vina-fail by category
cats = {"Biological TM\n(Zn Mg Mn Fe Co Ni Cu)": ["zn","mg","mn","fe","co","ni","cu"],
        "Medicinal\n(Pt Pd Ru Ir Au Rh Ag)": ["pt","pd","ru","ir","au","rh","ag"],
        "Toxicology\n(Cd Hg Tl Pb As Sb Bi)": ["cd","hg","tl","pb","as","sb","bi"],
        "s-block / early TM": ["na","k","li","al","sr","ba","v","cr","ti","sc","y","zr","nb","hf","ta","w","mo"],
        "Lanthanides /\nActinides": [r["mode"] for r in metal if r["mode"] in
            ("la","ce","pr","nd","sm","eu","gd","tb","dy","ho","er","tm","yb","lu",
             "ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","uo2","pm","ra")],
        "Other\n(metalloids etc.)": [r["mode"] for r in metal if r["mode"] not in
            ("zn","mg","mn","fe","co","ni","cu","pt","pd","ru","ir","au","rh","ag",
             "cd","hg","tl","pb","as","sb","bi","na","k","li","al","sr","ba","v","cr","ti","sc","y","zr","nb","hf","ta","w","mo",
             "la","ce","pr","nd","sm","eu","gd","tb","dy","ho","er","tm","yb","lu",
             "ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","uo2","pm","ra")]}
labels = list(cats.keys())
v_ok, v_fail = [], []
for lab, modes in cats.items():
    ok = sum(1 for m in modes for r in metal if r["mode"] == m and r["vina_rc"] == 0)
    fail = sum(1 for m in modes for r in metal if r["mode"] == m and r["vina_rc"] != 0)
    v_ok.append(ok); v_fail.append(fail)
y = np.arange(len(labels))
ax2.barh(y, v_ok, 0.4, color=C_LKINA, edgecolor="#333", linewidth=0.5, label="Vina 1.2.7 OK")
ax2.barh(y, v_fail, 0.4, left=v_ok, color="#E8D9D7", edgecolor="#333", linewidth=0.5, label="Vina rejects")
ax2.set_yticks(y); ax2.set_yticklabels(labels, fontsize=9)
ax2.set_xlabel("modes")
ax2.set_xlim(0, 40)
ax2.set_title("(b) Vina acceptance by metal family", fontsize=11.5)
ax2.spines[["top", "right"]].set_visible(False)
ax2.legend(loc="lower right", frameon=False, fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig5_metal_coverage_110.png"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 6 — coordination accuracy (left) + pseudoatom geometry (right)
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.3), gridspec_kw={"width_ratios": [1.2, 1]})

errs = [abs(r["d_lig"] - r["d0"]) for r in metal if r["d_lig"]]
bins = np.arange(0, 1.01, 0.1)
ax1.hist(errs, bins=bins, color=C_LKINA, edgecolor="white", alpha=0.92)
ax1.axvline(0.2, color=C_GOLD, linestyle="--", lw=1.2)
ax1.text(0.22, max(np.histogram(errs, bins=bins)[0]) * 0.9, "mean = 0.20 Å", color=C_GOLD, fontsize=10, rotation=90, va="top")
ax1.set_xlabel("|d(ligand–metal) − r_eq|  (Å)")
ax1.set_ylabel("metal modes")
ax1.set_title("(a) Coordination accuracy, 110 modes\n108/110 within 0.5 Å of r_eq", fontsize=11)
ax1.spines[["top", "right"]].set_visible(False)
ax1.grid(axis="y", color="#E5E7EB", linewidth=0.7, zorder=0)

# pseudoatom geometry check: measured vs expected donor distances per geometry
geo = ff["pseudoatom_geometry"]
labels = [f"{g['mode']} ({g['pseudo']})" for g in geo]
meas = [np.mean(g["dists"]) if g["dists"] else 0 for g in geo]
d0   = [g["d0"] for g in geo]
xx = np.arange(len(labels))
ax2.bar(xx - 0.18, meas, 0.36, color=C_LKINA, label="measured d")
ax2.bar(xx + 0.18, d0, 0.36, color=C_GREY, label="expected r_eq")
for i, g in enumerate(geo):
    ax2.text(i, max(meas[i], d0[i]) + 0.05, f"{g['n_ideal']}/{g['n_donors_checked']}", ha="center", fontsize=9, color="#333")
ax2.set_xticks(xx); ax2.set_xticklabels(labels, fontsize=9)
ax2.set_ylabel("donor–metal distance (Å)")
ax2.set_ylim(0, 3.0)
ax2.set_title("(b) Pseudoatom geometry check (TZ/SQ/MH/JT)\nideal-count / donors-checked", fontsize=11)
ax2.spines[["top", "right"]].set_visible(False)
ax2.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig6_coordination_geometry.png"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 7 — BVS oxidation-state inference (14/14)
# ---------------------------------------------------------------------------
bvs = ff["bvs_inference"]
labels = [f"{b['tok']}@{b['d0']:.2f} Å" for b in bvs]
colors = [C_LKINA if b["ok"] else "#C0392B" for b in bvs]
fig, ax = plt.subplots(figsize=(10.4, 3.8))
y = np.arange(len(bvs))
for i, b in enumerate(bvs):
    ax.barh(i, 1.0, 0.62, color=(C_LKINA if b["ok"] else "#E8D9D7"), edgecolor="#333", linewidth=0.5)
    ax.text(0.5, i, f"{b['expected']}", ha="center", va="center", fontsize=10,
            color="white" if b["ok"] else "#333", fontweight="bold")
    ax.text(1.03, i, f"→ {b['detected']}", va="center", fontsize=10,
            color=C_LKINA if b["ok"] else "#C0392B", fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5)
ax.set_xlim(0, 1.55)
ax.set_xlabel("Bond-Valence-Sum inference (receptor donors, auto-detected mode)")
ax.set_title("BVS oxidation-state inference: 14/14 correct  (Fe/Cu/Mn/Co/V/Mo/Ni, ±1 e⁻)", fontsize=11.5)
ax.spines[["top", "right"]].set_visible(False)
ax.set_xticks([])
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig7_bvs_inference.png"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 8 — covalent framework
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(10.4, 5.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.9, 1.35])

# (a) P1 / P1+P2 / P4 / C3 success per preset
ax1 = fig.add_subplot(gs[0])
presets = [p["preset"] for p in cov["presets"]]
rows = cov["presets"]
names = ["P1\ndistance", "P1+P2\nangle", "P4\nvdW-scale", "C3\ntwo-step"]
vals = [[1 if r["p1"]["rc"] == 0 else 0 for r in rows],
        [1 if r["p12"]["rc"] == 0 else 0 for r in rows],
        [1 if all(r["p4_vdw_scale"][k]["rc"] == 0 for k in ("0.0","0.5","1.0")) else 0 for r in rows],
        [1 if r["c3_two_step"]["rc"] == 0 else 0 for r in rows]]
x = np.arange(4)
ax1.bar(x, [sum(v) for v in vals], 0.55, color=[C_LKINA, C_GOLD, C_BLUE, "#7A9E7E"], edgecolor="#333", linewidth=0.6)
for xi, v in zip(x, [sum(v) for v in vals]):
    ax1.text(xi, v + 0.1, f"{v}/6", ha="center", fontweight="bold")
ax1.set_xticks(x); ax1.set_xticklabels(names, fontsize=9)
ax1.set_ylim(0, 7); ax1.set_ylabel("presets completed")
ax1.set_title("(a) Covalent tiers P1–P4 + C3\n(6 reaction presets)", fontsize=11)
ax1.spines[["top", "right"]].set_visible(False)
ax1.grid(axis="y", color="#E5E7EB", linewidth=0.7, zorder=0)

# (b) NAC discrimination
ax2 = fig.add_subplot(gs[1])
nac = [r["p12"].get("nac") for r in rows]
nac_y = [1 if n == "YES" else 0 for n in nac]
preset_short = [r["preset"].replace("_", "\n") for r in rows]
colors2 = [C_LKINA if n == "YES" else C_GREY for n in nac]
ax2.bar(np.arange(len(nac)), nac_y, 0.6, color=colors2, edgecolor="#333", linewidth=0.5)
for i, (n, yv) in enumerate(zip(nac, nac_y)):
    ax2.text(i, yv + 0.05, n, ha="center", fontweight="bold", fontsize=9)
ax2.set_xticks(np.arange(len(nac))); ax2.set_xticklabels(preset_short, fontsize=8)
ax2.set_ylim(0, 1.25); ax2.set_yticks([0, 1]); ax2.set_yticklabels(["NO", "YES"])
ax2.set_title("(b) NAC detection\n(geometric discrimination)", fontsize=11)
ax2.spines[["top", "right"]].set_visible(False)

# (c) reactive distance energy well (Gaussian)
ax3 = fig.add_subplot(gs[2])
scan = cov["well_scan"]["points"]
ds = [p["dist"] for p in scan]
es = [p["e_reactive"] for p in scan]
ax3.plot(ds, es, "o-", color=C_LKINA, lw=1.8, markersize=5)
ax3.axvline(1.82, color=C_GOLD, linestyle="--", lw=1.1)
ax3.text(1.84, -6.3, "r₀ = 1.82 Å (Cys SG)", color=C_GOLD, fontsize=9)
ax3.set_xlabel("ligand C – receptor SG distance (Å)")
ax3.set_ylabel("reactive distance energy (kcal/mol)")
ax3.set_title("(c) P1 distance potential, cys_michael\nGaussian well, depth = 10 kcal/mol", fontsize=11)
ax3.spines[["top", "right"]].set_visible(False)
ax3.grid(color="#E5E7EB", linewidth=0.7, zorder=0)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig8_covalent_framework.png"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 9 — inline maps + metal-as-ligand
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.1), gridspec_kw={"width_ratios": [1.1, 1]})

gm = cov["generate_maps"]
maps = gm["maps"]
ax1.barh([0], [gm["n_maps"]], 0.5, color=C_LKINA, edgecolor="#333")
ax1.text(gm["n_maps"] + 3, 0, f"{gm['n_maps']} maps", va="center", fontweight="bold", color=C_LKINA)
ax1.set_yticks([]); ax1.set_xlim(0, 130)
ax1.set_xlabel("affinity maps written (probe types)")
ax1.set_title("(a) Inline AG4 grid generator\n(--generate_maps, no autogrid4)", fontsize=11)
ax1.spines[["top", "right"]].set_visible(False)

lm = ff["metal_as_ligand"]
lab = [l["label"].replace("_", " ") for l in lm]
pen = [float(l["ligand_metal_geom"]) if l["ligand_metal_geom"] else 0 for l in lm]
ax2.bar(np.arange(2), pen, width=0.45, color=[C_LKINA, C_GOLD], edgecolor="#333")
for i, p in enumerate(pen):
    ax2.text(i, p + 0.03, f"{p:.3f}", ha="center", fontweight="bold")
ax2.set_xticks(np.arange(2)); ax2.set_xticklabels(lab)
ax2.set_ylabel("geometry penalty (kcal/mol)")
ax2.set_title("(b) Metal-as-ligand geometry QC\nPt(II) square-planar penalty", fontsize=11)
ax2.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig9_maps_metal_ligand.png"), bbox_inches="tight")
plt.close(fig)

print("Figures written:")
for f in ["fig5_metal_coverage_110.png", "fig6_coordination_geometry.png",
          "fig7_bvs_inference.png", "fig8_covalent_framework.png",
          "fig9_maps_metal_ligand.png"]:
    p = os.path.join(OUT, f)
    print(f"  {f}  {os.path.getsize(p)} bytes")
