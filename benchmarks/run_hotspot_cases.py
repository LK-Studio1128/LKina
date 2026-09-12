#!/usr/bin/env python3
"""Hotspot-case covalent/precision docking tests.

1. 6OIM (KRAS G12C + MOV acrylamide inhibitor):
   reactive docking with cys_michael preset, reactive_rec_atom = A:12:SG,
   reactive_lig_atom = C17 (Michael-acceptor beta carbon).
   Compare: vina scoring, ad4 non-reactive, ad4 reactive (distance/hybrid),
   two-step C3.

2. 3HS4 (CA II + AZM sulfonamide, catalytic Zn):
   zn metal mode + metal_bias to test zinc-binder coordination recovery.

Output: hotspot_results.json
"""
import json
import math
import subprocess
import time
from pathlib import Path

BASE = Path(__file__).parent
HOT = BASE / "hotspot_cases"
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
OUT = BASE / "hotspot_results.json"
CASE_DIR = HOT

M6 = json.load(open(HOT / "6OIM_meta.json"))
M3 = json.load(open(HOT / "3HS4_meta.json"))
MC3 = tuple(M3["metal_coords"][0])

# crystal reference: CYS12 SG in 6OIM, AZM S atom near Zn in 3HS4
SG12 = (-6.344, -3.260, 0.409)


def box(m, pad=8):
    c, s = m["center"], m["span"]
    return ["--center_x", f"{c[0]:.2f}", "--center_y", f"{c[1]:.2f}", "--center_z", f"{c[2]:.2f}",
            "--size_x", f"{max(s[0]+pad,12):.0f}", "--size_y", f"{max(s[1]+pad,12):.0f}",
            "--size_z", f"{max(s[2]+pad,12):.0f}"]


def dock(args, tag, timeout=1800):
    outp = CASE_DIR / f"{tag}.pdbqt"
    cmd = [LKINA, "--out", str(outp), "--verbosity", "1"] + args
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        dt = round(time.time() - t0, 2)
    except subprocess.TimeoutExpired:
        return {"tag": tag, "fail": "timeout"}
    res = {"tag": tag, "runtime_s": dt}
    es = [float(ln.split()[3]) for ln in open(outp, errors="ignore")
          if "VINA RESULT" in ln] if outp.exists() else []
    if es:
        res["best"] = es[0]
        res["energies"] = es[:5]
    else:
        res["fail"] = (r.stderr or "no output").strip()[-200:]
    return res


def first_model_coords(p):
    coords, in_m = [], -1
    try:
        for ln in open(p, errors="ignore"):
            if ln.startswith("MODEL"):
                in_m += 1
            elif ln.startswith(("ATOM", "HETATM")) and in_m == 0:
                coords.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            elif ln.startswith("ENDMDL") and in_m == 0:
                break
    except FileNotFoundError:
        return []
    return coords


def nearest_atom_dist(p, ref, names=None):
    """Min dist from docked model-0 atoms (optionally filtered by atom name) to ref point."""
    best = None
    in_m = -1
    for ln in open(p, errors="ignore"):
        if ln.startswith("MODEL"):
            in_m += 1
        elif ln.startswith(("ATOM", "HETATM")) and in_m == 0:
            nm = ln[12:16].strip()
            if names and nm not in names:
                continue
            c = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            d = math.dist(c, ref)
            if best is None or d < best:
                best = d
        elif ln.startswith("ENDMDL") and in_m == 0:
            break
    return round(best, 2) if best is not None else None


results = {"6OIM_covalent": [], "3HS4_zinc": []}

# ---------------- 6OIM covalent battery
REC6, LIG6 = str(HOT / "6OIM_rec.pdbqt"), str(HOT / "6OIM_lig.pdbqt")
BOX6 = box(M6)
jobs6 = [
    ("6OIM_vina",       ["--receptor", REC6, "--ligand", LIG6, "--scoring", "vina"] + BOX6 +
                        ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6"]),
    ("6OIM_ad4_std",    ["--receptor", REC6, "--ligand", LIG6, "--scoring", "ad4",
                         "--generate_maps", "--no_auto_metal"] + BOX6 +
                        ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6"]),
    ("6OIM_react_dist", ["--receptor", REC6, "--ligand", LIG6, "--scoring", "ad4",
                         "--generate_maps"] + BOX6 +
                        ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6",
                         "--reactive_mode", "distance", "--reactive_preset", "cys_michael",
                         "--reactive_rec_atom", "A:12:SG", "--reactive_lig_atom", "name:C17"]),
    ("6OIM_react_hyb",  ["--receptor", REC6, "--ligand", LIG6, "--scoring", "ad4",
                         "--generate_maps"] + BOX6 +
                        ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6",
                         "--reactive_mode", "hybrid", "--reactive_hybrid_vdw_scale", "0.3",
                         "--reactive_preset", "cys_michael",
                         "--reactive_rec_atom", "A:12:SG", "--reactive_lig_atom", "name:C17"]),
    ("6OIM_react_2step", ["--receptor", REC6, "--ligand", LIG6, "--scoring", "ad4",
                          "--generate_maps"] + BOX6 +
                         ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6",
                          "--reactive_mode", "distance", "--reactive_preset", "cys_michael",
                          "--reactive_rec_atom", "A:12:SG", "--reactive_lig_atom", "name:C17",
                          "--reactive_two_step", "1", "--reactive_presample_dist", "8"]),
]
print("== 6OIM covalent ==")
for tag, args in jobs6:
    r = dock(args, tag)
    p = CASE_DIR / f"{tag}.pdbqt"
    # C17 (ligand warhead C) to CYS12 SG distance = surrogate covalent-bond metric
    r["warhead_to_SG_dist"] = nearest_atom_dist(p, SG12, names={"C17", "SG"})
    r["warhead_C17_to_SG"] = nearest_atom_dist(p, SG12, names={"C17"})
    results["6OIM_covalent"].append(r)
    print(tag, r.get("best"), r.get("warhead_C17_to_SG"), r.get("runtime_s"), r.get("fail", "")[:80])

# ---------------- 3HS4 zinc battery
REC3, LIG3 = str(HOT / "3HS4_rec.pdbqt"), str(HOT / "3HS4_lig.pdbqt")
BOX3 = box(M3)
jobs3 = [
    ("3HS4_vina",     ["--receptor", REC3, "--ligand", LIG3, "--scoring", "vina"] + BOX3 +
                      ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6"]),
    ("3HS4_ad4_std",  ["--receptor", REC3, "--ligand", LIG3, "--scoring", "ad4",
                       "--generate_maps", "--no_auto_metal"] + BOX3 +
                      ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6"]),
    ("3HS4_zn",       ["--receptor", REC3, "--ligand", LIG3, "--scoring", "ad4",
                       "--generate_maps", "--metal_mode", "zn"] + BOX3 +
                      ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6"]),
    ("3HS4_zn_bias",  ["--receptor", REC3, "--ligand", LIG3, "--scoring", "ad4",
                       "--generate_maps", "--metal_mode", "zn", "--metal_bias",
                       "--metal_bias_strength", "2",
                       # NOTE(v1.0.1): metal_bias needs explicit ligand atom spec
                       "--reactive_lig_atom", "index:1"] + BOX3 +
                      ["--seed", "42", "--exhaustiveness", "8", "--num_modes", "3", "--cpu", "6"]),
]
print("== 3HS4 zinc ==")
for tag, args in jobs3:
    r = dock(args, tag)
    p = CASE_DIR / f"{tag}.pdbqt"
    r["min_dist_to_Zn"] = nearest_atom_dist(p, MC3)          # any ligand heavy atom
    r["S_to_Zn_dist"] = nearest_atom_dist(p, MC3, names={"S"})  # AZM sulfonamide S
    results["3HS4_zinc"].append(r)
    print(tag, r.get("best"), r.get("min_dist_to_Zn"), r.get("S_to_Zn_dist"), r.get("fail", "")[:80])

json.dump(results, open(OUT, "w"), indent=1)
print("saved ->", OUT)
