#!/usr/bin/env python3
"""Parameter sweep battery for LKina.

Sweeps (each on a representative system):
  A. exhaustiveness 1/2/4/8/16/32 x 3 seeds  -> success-rate + convergence + runtime
     system: 3HS4 (CA II Zn, ad4+zn mode)
  B. seed variance: seeds 1..5, exhaustiveness 8 -> energy/rmsd spread
     systems: 3HS4 (zn), 6OIM (vina scoring)
  C. metal_soft_weight 0/0.1/0.2/0.3/0.5 -> metal-ligand geometry (min Zn-lig dist)
     system: 3HS4
  D. metal_bias_strength 0/1/2/4/8 (+width 1.5) -> min Zn-lig dist + energy
     system: 3HS4, ad4
  E. reactive_attractor_strength 2/4/8/12/16 -> covalent anchor distance
     system: reactive_tests/cys_michael (lig+rec)

Output: parameter_sweep_results.json
"""
import json
import math
import subprocess
import time
from pathlib import Path

BASE = Path(__file__).parent
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
HOT = BASE / "hotspot_cases"
OUTJ = BASE / "parameter_sweep_results.json"

REC3 = str(HOT / "3HS4_rec.pdbqt")
LIG3 = str(HOT / "3HS4_lig.pdbqt")
REC6 = str(HOT / "6OIM_rec.pdbqt")
LIG6 = str(HOT / "6OIM_lig.pdbqt")
M3 = json.load(open(HOT / "3HS4_meta.json"))
M6 = json.load(open(HOT / "6OIM_meta.json"))

def box(m, pad=8):
    c, s = m["center"], m["span"]
    return ["--center_x", f"{c[0]:.2f}", "--center_y", f"{c[1]:.2f}", "--center_z", f"{c[2]:.2f}",
            "--size_x", f"{max(s[0]+pad,12):.0f}", "--size_y", f"{max(s[1]+pad,12):.0f}",
            "--size_z", f"{max(s[2]+pad,12):.0f}"]

BOX3, BOX6 = box(M3), box(M6)
MC3 = tuple(M3["metal_coords"][0])

SWEEP_DIR = BASE / "parameter_sweeps"
SWEEP_DIR.mkdir(exist_ok=True)

def run(args, tag, timeout=1800):
    outp = SWEEP_DIR / f"{tag}.pdbqt"
    cmd = [LKINA, "--out", str(outp), "--verbosity", "0"] + args
    try:
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        dt = round(time.time() - t0, 2)
    except subprocess.TimeoutExpired:
        return {"tag": tag, "fail": "timeout"}
    res = {"tag": tag, "runtime_s": dt}
    es = []
    for ln in open(outp, errors="ignore") if outp.exists() else []:
        if "VINA RESULT" in ln:
            es.append(float(ln.split()[3]))
    if es:
        res["best"] = es[0]
        res["energies"] = es[:5]
    else:
        res["fail"] = (r.stderr or "no output").strip()[-150:]
    return res

def min_zn_dist(pdbqt, mc=MC3, model=0):
    coords = []
    in_model = -1
    try:
        for ln in open(pdbqt, errors="ignore"):
            if ln.startswith("MODEL"):
                in_model += 1
            elif ln.startswith(("ATOM", "HETATM")) and in_model == model:
                coords.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            elif ln.startswith("ENDMDL") and in_model == model:
                break
    except FileNotFoundError:
        return None
    if not coords:
        return None
    return round(min(math.dist(c, mc) for c in coords), 2)

def base3(extra):
    return ["--receptor", REC3, "--ligand", LIG3, "--scoring", "ad4",
            "--generate_maps"] + BOX3 + extra + ["--cpu", "6"]

def base6(extra):
    return ["--receptor", REC6, "--ligand", LIG6, "--scoring", "vina"] + BOX6 + \
           extra + ["--cpu", "6"]

results = {"A_exhaustiveness": [], "B_seed_variance": [], "C_soft_weight": [],
           "D_metal_bias": [], "E_reactive_strength": []}

# ---------------- A. exhaustiveness sweep (3HS4, zn mode, ad4)
print("== A exhaustiveness ==")
for exh in (1, 2, 4, 8, 16, 32):
    for seed in (42, 7, 2026):
        tag = f"A_exh{exh}_s{seed}"
        r = run(base3(["--metal_mode", "zn", "--seed", str(seed), "--exhaustiveness", str(exh),
                       "--num_modes", "1"]), tag)
        r["exhaustiveness"] = exh; r["seed"] = seed
        r["min_zn_dist"] = min_zn_dist(SWEEP_DIR / f"{tag}.pdbqt")
        results["A_exhaustiveness"].append(r)
        print(tag, r.get("best"), r.get("min_zn_dist"), r.get("runtime_s"), r.get("fail", ""))

# ---------------- B. seed variance (5 seeds), two systems
print("== B seed variance ==")
for seed in (1, 2, 3, 4, 5):
    tag = f"B_3HS4_s{seed}"
    r = run(base3(["--metal_mode", "zn", "--seed", str(seed), "--exhaustiveness", "8",
                   "--num_modes", "1"]), tag)
    r["system"] = "3HS4"; r["seed"] = seed
    r["min_zn_dist"] = min_zn_dist(SWEEP_DIR / f"{tag}.pdbqt")
    results["B_seed_variance"].append(r)
    print(tag, r.get("best"), r.get("min_zn_dist"))
for seed in (1, 2, 3, 4, 5):
    tag = f"B_6OIM_s{seed}"
    r = run(base6(["--seed", str(seed), "--exhaustiveness", "8", "--num_modes", "1"]), tag)
    r["system"] = "6OIM"; r["seed"] = seed
    results["B_seed_variance"].append(r)
    print(tag, r.get("best"))

# ---------------- C. metal_soft_weight sweep
print("== C soft_weight ==")
for w in (0.0, 0.1, 0.2, 0.3, 0.5):
    tag = f"C_sw{w}"
    r = run(base3(["--metal_mode", "zn", "--seed", "42", "--exhaustiveness", "8",
                   "--num_modes", "1", "--metal_soft_weight", str(w)]), tag)
    r["soft_weight"] = w
    r["min_zn_dist"] = min_zn_dist(SWEEP_DIR / f"{tag}.pdbqt")
    results["C_soft_weight"].append(r)
    print(tag, r.get("best"), r.get("min_zn_dist"))

# ---------------- D. metal_bias sweep (bias needs generate_maps; strength on attractor)
print("== D metal_bias ==")
for s in (0.0, 1.0, 2.0, 4.0, 8.0):
    tag = f"D_mb{s}"
    r = run(base3(["--metal_mode", "zn", "--seed", "42", "--exhaustiveness", "8",
                   "--num_modes", "1", "--metal_bias", "--metal_bias_strength", str(s),
                   # NOTE(v1.0.1): O5 metal_bias injects rec anchor but does not
                   # auto-fill the ligand atom spec; pass index:1 to activate.
                   "--reactive_lig_atom", "index:1"]), tag)
    r["bias_strength"] = s
    r["min_zn_dist"] = min_zn_dist(SWEEP_DIR / f"{tag}.pdbqt")
    results["D_metal_bias"].append(r)
    print(tag, r.get("best"), r.get("min_zn_dist"))

# ---------------- E. reactive attractor strength sweep (cys_michael synthetic)
print("== E reactive strength ==")
RT = BASE / "reactive_tests" / "cys_michael"
RECR, LIGR = str(RT / "rec.pdbqt"), str(RT / "lig.pdbqt")
for st in (2, 4, 8, 12, 16):
    tag = f"E_rs{st}"
    args = ["--receptor", RECR, "--ligand", LIGR, "--scoring", "ad4", "--generate_maps",
            "--center_x", "0", "--center_y", "0", "--center_z", "0",
            "--size_x", "20", "--size_y", "20", "--size_z", "20",
            "--seed", "42", "--exhaustiveness", "8", "--num_modes", "1", "--cpu", "6",
            "--reactive_mode", "distance", "--reactive_preset", "cys_michael",
            "--reactive_rec_atom", "A:1:SG", "--reactive_lig_atom", "index:1",
            "--reactive_attractor_strength", str(st)]
    r = run(args, tag)
    r["attractor_strength"] = st
    results["E_reactive_strength"].append(r)
    print(tag, r.get("best"), r.get("runtime_s"), r.get("fail", ""))

json.dump(results, open(OUTJ, "w"), indent=1)
print("saved ->", OUTJ)
