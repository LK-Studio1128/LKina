#!/usr/bin/env python3
"""Generate publication-quality figures for the LKina manuscript (v3).

Journal-compliance upgrades over v2:
  * 300-dpi PNG (print) + vector PDF (submission) for every figure
  * Arial font family, panel labels (a)/(b)/(c) as bold text outside axes
  * consistent muted colorblind-safe palette; no seaborn styling
  * axis units in parentheses; tick direction out; single-stroke spines
  * error-bar / annotation typography >= 8 pt at final size (single col 3.5in, double 7.2in)

New supplementary figures:
  Fig S1  metal-mode energy-well landscape across all 110 modes
          (E_ideal < E_far coordination-well confirmation)
  Fig S2  NAC attack-angle discrimination: measured angle vs preset target
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
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        9,
    "axes.edgecolor":  "#333333",
    "axes.linewidth":   0.8,
    "axes.titlesize":  10,
    "axes.labelsize":   9.5,
    "legend.fontsize":  8,
    "xtick.labelsize":  8.5,
    "ytick.labelsize":  8.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi":      600,
    "savefig.dpi":     600,
})

# colorblind-safe-ish palette, print friendly
C_LKINA = "#B2182B"   # dark red
C_VINA  = "#8FA3B3"   # grey blue
C_GREY  = "#C9CDD1"
C_GOLD  = "#D18A00"
C_BLUE  = "#2166AC"
C_GREEN = "#5A8F5A"

def style_ax(ax):
    ax.grid(axis="y", color="#E3E6E8", linewidth=0.6, zorder=0)

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".tif"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name, "(png + pdf + tif)")

metal = json.load(open(f"{BENCH}/metal_coverage_results_all.json"))
ff    = json.load(open(f"{BENCH}/feature_family_results.json"))
cov   = json.load(open(f"{BENCH}/covalent_full_results.json"))
rp    = json.load(open(f"{BENCH}/reactive_presets_results.json"))

# ---------------------------------------------------------------------------
# Figure 5 — coverage 110 modes: LKina vs Vina (a) + by-family acceptance (b)
# ---------------------------------------------------------------------------
n_ok   = sum(1 for r in metal if r["lkina_rc"] == 0)
n_vok  = sum(1 for r in metal if r["vina_rc"] == 0)
n_fail = 110 - n_vok

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1, 1.35]})

x = np.arange(2)
ax1.bar(x, [n_ok, n_vok], 0.52, color=[C_LKINA, C_VINA], edgecolor="#222", linewidth=0.6, zorder=3)
ax1.text(0, n_ok + 3, f"{n_ok}/110", ha="center", fontweight="bold", color=C_LKINA, fontsize=11)
ax1.text(1, n_vok + 3, f"{n_vok}/110", ha="center", fontweight="bold", color="#555", fontsize=11)
ax1.annotate("76 metal types rejected at PDBQT parse:\n\u201cAtom type X is not a valid AutoDock type\u201d",
             xy=(1.05, 40), xytext=(0.44, 135), fontsize=7.5, color="#333", va="top",
             arrowprops=dict(arrowstyle="->", color="#666", lw=0.8))
ax1.set_xticks(x); ax1.set_xticklabels(["LKina", "Vina 1.2.7"])
ax1.set_ylabel("Metal modes docked successfully")
ax1.set_ylim(0, 145)
ax1.set_title("Coverage of the complete metal_mode enum (110)", fontsize=9.5)
style_ax(ax1)
ax1.text(-0.16, 1.16, "A", transform=ax1.transAxes, fontsize=13, fontweight="bold", clip_on=False)

fam = {
    "Biological TM":        ("zn mg mn fe co ni cu".split(), None),
    "Medicinal":            ("pt pd ru ir au rh ag".split(), None),
    "Toxicology":           ("cd hg tl pb as sb bi".split(), None),
    "s-block / p-block":    ("na k li al sr ba ga in sn be b se ge te po at".split(), None),
    "Early TM":             ("v cr ti sc y zr nb hf ta w mo".split(), None),
    "Lanthanide /\nactinide": (["la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb","lu",
                                 "ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","ra","uo2"], None),
}
labels, vok_l, vfail_l = [], [], []
for lab, (modes, _) in fam.items():
    ok = sum(1 for m in modes for r in metal if r["mode"] == m and r["vina_rc"] == 0)
    fail = sum(1 for m in modes for r in metal if r["mode"] == m and r["vina_rc"] != 0)
    labels.append(lab); vok_l.append(ok); vfail_l.append(fail)
y = np.arange(len(labels))
ax2.barh(y, vok_l, 0.58, color=C_LKINA, edgecolor="#222", linewidth=0.5, label="Vina completes", zorder=3)
ax2.barh(y, vfail_l, 0.58, left=vok_l, color="#EBDBD9", edgecolor="#222", linewidth=0.5,
         label="Vina rejects (parse)", zorder=3)
for i, (o, f_) in enumerate(zip(vok_l, vfail_l)):
    ax2.text(o + f_ + 0.4, i, f"{o}/{o+f_}", va="center", fontsize=7.5, color="#333")
ax2.set_yticks(y); ax2.set_yticklabels(labels, fontsize=8)
ax2.set_xlabel("Number of metal modes")
ax2.set_xlim(0, 36)
ax2.invert_yaxis()
ax2.legend(loc="upper right", frameon=False, fontsize=7, borderaxespad=0.1)
ax2.grid(axis="x", color="#E3E6E8", linewidth=0.6, zorder=0)
ax2.text(-0.14, 1.16, "B", transform=ax2.transAxes, fontsize=13, fontweight="bold", clip_on=False)

fig.tight_layout(w_pad=2.4)
save(fig, "fig1_metal_coverage_110")

# ---------------------------------------------------------------------------
# Figure 6 — coordination accuracy histogram (a) + geometry check (b)
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.15, 1]})

errs = [abs(r["d_lig"] - r["d0"]) for r in metal if r["d_lig"]]
bins = np.arange(0, 1.01, 0.1)
cnt, edges, patches = ax1.hist(errs, bins=bins, color=C_LKINA, edgecolor="white", linewidth=0.7, zorder=3)
ax1.axvline(0.20, color=C_GOLD, linestyle="--", lw=1.1, zorder=4)
ax1.text(0.215, max(cnt)*0.40, "mean = 0.20 \u00C5", color=C_GOLD, fontsize=8.5, rotation=90, va="bottom",
         bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
for e, p in zip(edges[:-1], patches):
    if e < 0.5:
        p.set_alpha(1.0)
ax1.axvspan(0, 0.5, color=C_BLUE, alpha=0.06, zorder=1)
ax1.text(0.25, max(cnt)*1.10, "108/110 within 0.5 \u00C5", ha="center", fontsize=8, color=C_BLUE)
ax1.set_ylim(0, max(cnt)*1.24)
ax1.set_xlabel("|d(ligand\u2013metal) \u2212 $d_0$| (\u00C5)")
ax1.set_ylabel("Metal modes")
ax1.set_xlim(0, 1.0)
ax1.set_title("Coordination-distance accuracy (110 modes)", fontsize=9.5)
style_ax(ax1)
ax1.text(-0.15, 1.16, "A", transform=ax1.transAxes, fontsize=13, fontweight="bold", clip_on=False)

geo = ff["pseudoatom_geometry"]
labels_g = [f"{g['mode']}\n({g['pseudo']})" for g in geo]
meas = [np.mean(g["dists"]) for g in geo]
d0s  = [g["d0"] for g in geo]
xx = np.arange(len(labels_g))
b1 = ax2.bar(xx - 0.19, meas, 0.38, color=C_LKINA, label="measured d(M\u2013donor)", zorder=3)
b2 = ax2.bar(xx + 0.19, d0s,  0.38, color=C_GREY, label="expected r_eq", zorder=3)
for i, g in enumerate(geo):
    ax2.text(i, max(meas[i], d0s[i]) + 0.07, f"{g['n_ideal']}/{g['n_donors_checked']}",
             ha="center", fontsize=8, color="#333")
    ax2.errorbar(i - 0.19, meas[i], yerr=max(np.std(g["dists"]), 0.01), fmt="none",
                 ecolor="#222", elinewidth=0.8, capsize=2.5, zorder=4)
ax2.set_xticks(xx); ax2.set_xticklabels(labels_g, fontsize=8)
ax2.set_ylabel("Donor\u2013metal distance (\u00C5)")
ax2.set_ylim(0, 3.1)
ax2.legend(frameon=False, fontsize=7, loc="upper center", ncol=2, borderaxespad=0.1)
ax2.set_title("--metal_geometry_check (TZ/SQ/MH/JT)", fontsize=9.5)
style_ax(ax2)
ax2.text(-0.17, 1.16, "B", transform=ax2.transAxes, fontsize=13, fontweight="bold", clip_on=False)

fig.tight_layout(w_pad=2.2)
save(fig, "fig2_coordination_geometry")

# ---------------------------------------------------------------------------
# Figure 7 — BVS inference (14/14)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(3.5, 3.6))
bvs = ff["bvs_inference"]
labels_b = [f"{b['tok']} {b['d0']:.2f}\u00C5" for b in bvs]
y = np.arange(len(bvs))
for i, b in enumerate(bvs):
    ax.barh(i, 1.0, 0.62, color=(C_LKINA if b["ok"] else C_GREY), edgecolor="#222", linewidth=0.5, zorder=3)
    ax.text(0.5, i, b["expected"], ha="center", va="center", fontsize=8,
            color=("white" if b["ok"] else "#333"), fontweight="bold")
    ax.text(1.05, i, f"\u2192 {b['detected']}", va="center", fontsize=8,
            color=(C_LKINA if b["ok"] else "#B2182B"), fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels(labels_b, fontsize=8)
ax.invert_yaxis()
ax.set_xlim(0, 1.6)
ax.set_xlabel("Designed state \u2192 auto-detected state (BVS)")
ax.set_xticks([])
ax.set_title("BVS oxidation-state inference:\n14/14 correct (Fe/Cu/Mn/Co/V/Mo/Ni, $\\pm 1\\,e^-$)", fontsize=9.5)
fig.tight_layout()
save(fig, "fig3_bvs_inference")

# ---------------------------------------------------------------------------
# Figure 8 — covalent framework: tiers (a) NAC (b) well scan (c)
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 3.1))
gs = fig.add_gridspec(1, 3, width_ratios=[0.95, 1.25, 1.15])

ax1 = fig.add_subplot(gs[0])
rows = cov["presets"]
names = ["P1\ndist.", "P1+P2\nangle", "P4\nvdW", "C3\ntwo-step"]
vals = [[sum(1 for r in rows if r["p1"]["rc"] == 0)],
        ]
per_tier = [
    sum(1 for r in rows if r["p1"]["rc"] == 0),
    sum(1 for r in rows if r["p12"]["rc"] == 0),
    sum(1 for r in rows if all(r["p4_vdw_scale"][k]["rc"] == 0 for k in ("0.0","0.5","1.0"))),
    sum(1 for r in rows if r["c3_two_step"]["rc"] == 0),
]
x = np.arange(4)
bars = ax1.bar(x, per_tier, 0.55, color=[C_LKINA, C_GOLD, C_BLUE, "#4A7A4A"], edgecolor="#222",
               linewidth=0.6, zorder=3)
for xi, v in zip(x, per_tier):
    ax1.text(xi, v + 0.12, f"{v}/6", ha="center", fontweight="bold", fontsize=9)
ax1.set_xticks(x); ax1.set_xticklabels(names, fontsize=7.5)
ax1.set_ylim(0, 7.2); ax1.set_ylabel("Presets completed")
ax1.set_title("Four-tier framework, six presets", fontsize=9.5)
ax1.grid(axis="y", color="#E3E6E8", linewidth=0.6, zorder=0)
ax1.text(-0.22, 1.16, "A", transform=ax1.transAxes, fontsize=13, fontweight="bold", clip_on=False)

# (b) measured P1+P2 attack angle vs preset target
ax2 = fig.add_subplot(gs[1])
th0_map = {"cys_michael":109.5, "cys_sn2":180.0, "ser_covalent":109.5,
           "lys_targeting":109.5, "boronic_acid":None, "tyr_covalent":109.5}
short = {"cys_michael":"Michael\n(Cys)", "cys_sn2":"SN2\n(Cys)", "ser_covalent":"acylation\n(Ser)",
         "lys_targeting":"Schiff\n(Lys)", "boronic_acid":"boronate\n(Ser/OH)", "tyr_covalent":"attack\n(Tyr)"}
yy = np.arange(len(rows))
for i, r in enumerate(rows):
    t0 = th0_map[r["preset"]]
    ang = r["p12"]["angle"]
    ok = (r["p12"]["nac"] == "YES")
    ax2.plot(ang, yy[i], "o", ms=7, color=(C_LKINA if ok else C_GREY),
             markeredgecolor="#222", markeredgewidth=0.5, zorder=4)
    ax2.text(ang, yy[i] - 0.34, f"{ang:.0f}\u00B0", ha="center", fontsize=6.8, color="#444")
    if t0 is not None:
        ax2.plot([t0], yy[i], marker="*", ms=11, color=C_GOLD, markeredgecolor="#222",
                 markeredgewidth=0.4, zorder=4)
        half = 25.0 if t0 != 180.0 else 15.0
        ax2.barh(yy[i], 2*half, left=t0-half, height=0.5, color=C_GOLD, alpha=0.12, zorder=2)
ax2.axvline(np.nan, color=C_GOLD)
ax2.plot([], [], "*", ms=11, color=C_GOLD, markeredgecolor="#222",
         label="preset target $\\theta_0$ (\u00B125\u00B0)")
ax2.plot([], [], "o", ms=7, color=C_LKINA, markeredgecolor="#222", label="measured P1+P2 angle")
ax2.plot([], [], "o", ms=7, color=C_GREY, markeredgecolor="#222", label="NAC = NO")
ax2.set_yticks(yy); ax2.set_yticklabels([short[r["preset"]] for r in rows], fontsize=7.5)
ax2.invert_yaxis()
ax2.set_xlim(70, 205)
ax2.set_xlabel(r"P1+P2 attack angle ($^\circ$)")
ax2.legend(frameon=False, fontsize=6.6, loc="upper center", bbox_to_anchor=(0.5, 1.42), ncol=1)
ax2.set_title("NAC geometric discrimination", fontsize=9.5)
ax2.grid(axis="x", color="#E3E6E8", linewidth=0.6, zorder=0)
ax2.text(-0.28, 1.16, "B", transform=ax2.transAxes, fontsize=13, fontweight="bold", clip_on=False)

# (c) Gaussian distance well
ax3 = fig.add_subplot(gs[2])
scan = cov["well_scan"]["points"]
ds = [p["dist"] for p in scan]; es = [p["e_reactive"] for p in scan]
tot = [p["e_total"] for p in scan]
ax3.plot(ds, es, "o-", color=C_LKINA, lw=1.6, markersize=4.5, zorder=3, label="reactive term")
ax3.plot(ds, tot, "s--", color=C_VINA, lw=1.1, markersize=3.5, zorder=3, label="total energy")
ax3.axvline(1.82, color=C_GOLD, linestyle="--", lw=1.0, zorder=2)
ax3.annotate(r"$r_0$ = 1.82 Å", xy=(1.82, 500), xytext=(1.98, 620),
             fontsize=8, color=C_GOLD, arrowprops=dict(arrowstyle="->", color=C_GOLD, lw=0.8))
ax3.set_xlabel("Ligand C \u2013 receptor SG distance (\u00C5)")
ax3.set_ylabel("Energy (kcal/mol)")
ax3.legend(frameon=False, fontsize=7)
ax3.set_title("P1 distance potential (cys_michael)\nGaussian well depth 10 kcal/mol", fontsize=9.5)
ax3.grid(color="#E3E6E8", linewidth=0.6, zorder=0)
ax3.text(-0.24, 1.16, "C", transform=ax3.transAxes, fontsize=13, fontweight="bold", clip_on=False)

fig.tight_layout(w_pad=2.6)
save(fig, "fig5_covalent_framework")

# ---------------------------------------------------------------------------
# Figure 9 — maps audit (a) + metal-as-ligand QC (b)
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.1, 1]})

gm = cov["generate_maps"]
ax1.barh([0], [gm["n_maps"]], 0.42, color=C_LKINA, edgecolor="#222", zorder=3)
ax1.text(gm["n_maps"] + 2, 0, f"{gm['n_maps']} .map files", va="center",
         fontweight="bold", color=C_LKINA, fontsize=10)
ax1.set_yticks([]); ax1.set_xlim(0, 132)
ax1.set_xlabel("AutoGrid4-format affinity maps written in one pass")
ax1.set_title("--generate_maps (no external autogrid4)", fontsize=9.5)
ax1.grid(axis="x", color="#E3E6E8", linewidth=0.6, zorder=0)
ax1.text(-0.18, 1.16, "A", transform=ax1.transAxes, fontsize=13, fontweight="bold", clip_on=False)

lm = ff["metal_as_ligand"]
lab = [l["label"].replace("_", "\n") for l in lm]
pen = [float(l["ligand_metal_geom"]) for l in lm]
ax2.bar(np.arange(2), pen, width=0.42, color=[C_LKINA, C_GOLD], edgecolor="#222", zorder=3)
for i, p in enumerate(pen):
    ax2.text(i, p + 0.03, f"{p:.3f}", ha="center", fontweight="bold", fontsize=9)
ax2.annotate("\u00D7300", xy=(1, pen[1]*0.55), xytext=(0.32, pen[1]*0.75),
             fontsize=9.5, fontweight="bold", color="#333",
             arrowprops=dict(arrowstyle="->", lw=0.9, color="#333"))
ax2.set_xticks(np.arange(2)); ax2.set_xticklabels(["ideal 90\u00B0\nsquare planar", "distorted 60\u00B0"])
ax2.set_ylabel("Geometry penalty (kcal/mol)")
ax2.set_ylim(0, 1.05)
ax2.set_title("Metal-as-ligand geometry QC, Pt(II)", fontsize=9.5)
style_ax(ax2)
ax2.text(-0.2, 1.16, "B", transform=ax2.transAxes, fontsize=13, fontweight="bold", clip_on=False)

fig.tight_layout(w_pad=2.6)
save(fig, "fig4_maps_metal_ligand")

# ---------------------------------------------------------------------------
# Figure S1 — NEW: coordination-energy well landscape over 110 modes
# ---------------------------------------------------------------------------
wells = [(r["mode"], r["pseudo"], r["e_far"] - r["e_ideal"]) for r in metal
         if r.get("e_ideal") is not None and r.get("e_far") is not None]
fig, ax = plt.subplots(figsize=(7.2, 2.9))
order = {"TZ":0, "SQ":1, "LIN":2, "MH":3, "JT":4}
xs = np.arange(len(wells))
ps = [order[w[1]] for w in wells]
vals_w = np.array([w[2] for w in wells])
cols = [C_LKINA if v > 0 else C_GREY for v in vals_w]
ax.scatter(xs, ps, s=np.clip(vals_w, 2, None)*2.2 + 6, c=cols, alpha=0.85,
           edgecolors="#222", linewidths=0.35, zorder=3)
ax.scatter([], [], s=60, c=C_LKINA, edgecolors="#222", linewidths=0.4, label="coordination well confirmed (104/110)")
ax.scatter([], [], s=30, c=C_GREY, edgecolors="#222", linewidths=0.4, label="no well at this probe position (6/110)")
ax.axhline(0, color="#888", lw=0.7, linestyle=":")
ax.set_yticks(range(5)); ax.set_yticklabels(list(order.keys()))
ax.invert_yaxis()
ax.set_xlabel("Synthetic metal-coordination systems (110 modes, sorted by engine enum order)")
ax.set_ylabel("Pseudoatom class")
ax.legend(frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(0.02, 0.06))
ax.set_title("Coordination potential well $E_{\\mathrm{ideal}} < E_{\\mathrm{far}}$ across the complete mode enum", fontsize=9.5)
fig.tight_layout()
save(fig, "figS1_energy_well_landscape")

# ---------------------------------------------------------------------------
# Figure S2 — NEW: reactive-preset P1 vs P1+P2 convergence
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(3.5, 2.9))
names_s = [r["preset"].split("_")[0] + "\n" + r["preset"].split("_")[1].title() for r in rp]
xx = np.arange(len(rp))
p1d  = [r["p1_dist"] for r in rp]
p12d = [r["p12_dist"] for r in rp]
r0m  = [{"cys_michael":1.82,"cys_sn2":1.82,"ser_covalent":1.34,
         "lys_targeting":1.47,"boronic_acid":1.47,"tyr_covalent":1.38}[r["preset"]] for r in rp]
ax.plot(xx, p1d, "o--", color=C_VINA, ms=6, label="P1 (distance only)", zorder=3)
ax.plot(xx, p12d, "o-", color=C_LKINA, ms=6, label="P1+P2 (+angle)", zorder=4)
ax.plot([], [], "_", color=C_GOLD, lw=2.4, label=r"preset $r_0$")
for xi, rv in zip(xx, r0m):
    ax.hlines(rv, xi-0.24, xi+0.24, color=C_GOLD, lw=2.4, zorder=5)
ax.set_xticks(xx); ax.set_xticklabels(names_s, fontsize=7.2)
ax.set_ylabel("Reactive distance (\u00C5)")
ax.legend(frameon=False, fontsize=7)
ax.set_title("Adding the angular constraint (P2)\npulls poses toward the ideal bond length", fontsize=9)
fig.tight_layout()
save(fig, "figS2_reactive_convergence")

print("All figures done (600 dpi PNG + vector PDF).")
