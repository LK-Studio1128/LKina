#!/usr/bin/env python3
"""Scoring-function comparison: vina / vinardo / ad4 on prepared systems.

Systems (from existing benchmarks + hotspots):
  4JC  (metalloprotein-ligand, Zn, from real_systems)
  3HS4 (CA II + AZM, Zn catalytic)   -- hotspot
  6OIM (KRAS G12C + MOV)             -- hotspot

For each system x scoring fn: dock, record energy, time, top pose -> metal
distance metrics. Output: scoring_comparison_results.json
"""
import json
import subprocess
import time
from pathlib import Path

BASE = Path(__file__).parent
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
HOT = BASE / "hotspot_cases"
REAL = BASE / "real_systems"
OUT = BASE / "scoring_comparison"

SEED = 42
CPU = 6

SYSTEMS = [
    # id, rec, lig, metal_mode, extra
    ("4JC",  REAL / "4JC_core_rec.pdbqt", BASE / "4JC305_lig.pdbqt", "zn", []),
    ("3HS4", HOT / "3HS4_rec.pdbqt", HOT / "3HS4_lig.pdbqt", "zn", []),
    ("6OIM", HOT / "6OIM_rec.pdbqt", HOT / "6OIM_lig.pdbqt", None, []),  # no metal in pocket
]

def meta_box(pdb_id):
    if pdb_id == "4JC":
        # 4JC box: pocket center used in run_feature_tests.py
        return ["--center_x", "-2.55", "--center_y", "2.29", "--center_z", "85.60",
                "--size_x", "25", "--size_y", "25", "--size_z", "25"]
    m = json.load(open(HOT / f"{pdb_id}_meta.json"))
    c, s = m["center"], m["span"]
    return ["--center_x", f"{c[0]:.2f}", "--center_y", f"{c[1]:.2f}", "--center_z", f"{c[2]:.2f}",
            "--size_x", f"{max(s[0]+8,12):.0f}", "--size_y", f"{max(s[1]+8,12):.0f}",
            "--size_z", f"{max(s[2]+8,12):.0f}"]

def read_energies(pdbqt):
    out = []
    for ln in open(pdbqt, errors="ignore"):
        if "REMARK VINA RESULT" in ln:
            parts = ln.split()
            out.append(float(parts[3]))
    return out

def read_coords(pdbqt, model=0):
    """Return heavy-atom coords of model N."""
    models, cur = [], []
    in_model = -1
    for ln in open(pdbqt, errors="ignore"):
        if ln.startswith("MODEL"):
            in_model += 1
            cur = []
        elif ln.startswith(("ATOM", "HETATM")):
            if in_model <= model:
                cur.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
        elif ln.startswith("ENDMDL"):
            if in_model == model:
                models.append(cur)
                break
    return models[0] if models else []

def metal_center(pdb_id):
    if pdb_id == "4JC":
        # ZN from core receptor (ZN A 302, per run_feature_tests.py setup)
        for ln in open(REAL / "4JC_core_rec.pdbqt", errors="ignore"):
            if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "ZN":
                return (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        return None
    m = json.load(open(HOT / f"{pdb_id}_meta.json"))
    return tuple(m["metal_coords"][0])

def closest_metal_dist(coords, mc):
    if mc is None or not coords:
        return None
    return min(((x-mc[0])**2 + (y-mc[1])**2 + (z-mc[2])**2) ** 0.5 for x, y, z in coords)

def main():
    OUT.mkdir(exist_ok=True)
    results = []
    for pdb_id, rec, lig, metal_mode, extra in SYSTEMS:
        box = meta_box(pdb_id)
        mc = metal_center(pdb_id)
        for scoring in ("vina", "vinardo", "ad4"):
            tag = f"{pdb_id}_{scoring}"
            outp = OUT / f"{tag}_out.pdbqt"
            maps_pref = OUT / f"maps_{tag}"
            # vina/vinardo: receptor+ligand mode (maps auto-computed)
            # ad4: needs --generate_maps from receptor+ligand
            cmd = [LKINA, "--receptor", str(rec), "--ligand", str(lig)] + box + [
                "--scoring", scoring,
                "--out", str(outp), "--seed", str(SEED), "--exhaustiveness", "8",
                "--num_modes", "5", "--cpu", str(CPU), "--verbosity", "1"]
            if scoring == "ad4":
                r = subprocess.run(cmd + ["--generate_maps", "--metal_mode", "zn"]
                                   if metal_mode else cmd + ["--generate_maps"],
                                   capture_output=True, text=True, timeout=3600)
            else:
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            t0 = time.time()
            # energy from output
            entry = {"id": pdb_id, "scoring": scoring, "metal_mode": metal_mode or "-"}
            if outp.exists():
                es = read_energies(outp)
                entry["energies"] = es
                entry["best"] = es[0] if es else None
                coords = read_coords(outp)
                entry["min_metal_dist_A"] = closest_metal_dist(coords, mc)
                entry["n_models"] = len(es)
            else:
                entry["fail"] = (r.stderr or "").strip()[-200:]
            # timing run (separate, verbosity 0, num_modes 1 for clean timing)
            outp2 = OUT / f"{tag}_time.pdbqt"
            tcmd = [LKINA, "--receptor", str(rec), "--ligand", str(lig)] + box + [
                "--scoring", scoring,
                "--out", str(outp2), "--seed", str(SEED), "--exhaustiveness", "8",
                "--num_modes", "1", "--cpu", str(CPU), "--verbosity", "0"]
            if scoring == "ad4":
                tcmd += ["--generate_maps"] + (["--metal_mode", "zn"] if metal_mode else [])
            try:
                t1 = time.time()
                subprocess.run(tcmd, capture_output=True, text=True, timeout=3600)
                entry["runtime_s"] = round(time.time() - t1, 2)
            except subprocess.TimeoutExpired:
                entry["runtime_s"] = None
            results.append(entry)
            print(json.dumps(entry, ensure_ascii=False))
    json.dump(results, open(BASE / "scoring_comparison_results.json", "w"), indent=1)
    print("saved -> scoring_comparison_results.json")

if __name__ == "__main__":
    main()
