#!/usr/bin/env python3
"""metal_soft_weight ablation — search-time soft constraint vs post-hoc rerank.

Paper Section 4.6 future-work item: systematic calibration of
--metal_soft_weight across Zn/Fe/Cu systems. Three real metalloprotein
systems, fixed grid, w in {0, 0.1, 0.3, 0.5}:

  w=0    : metal geometry terms only in post-docking rerank (default, paper §3.6)
  w>0    : log-sum-exp-smoothed metal geometry enters the search gradient
           (--metal_soft_weight), so the search is steered toward
           coordination-competent poses.

Metric: donor-metal distance of the best pose (same as Section 3.6). The
hypothesis is w>0 improves the fraction of poses placing a donor within 3.0 A
of the metal without hurting the (already weak) global RMSD.

Usage: python3 benchmarks/run_soft_weight_ablation.py
Output: benchmarks/soft_weight_results.json
"""
import json, os, re, subprocess, sys, glob
from pathlib import Path

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LKINA = os.path.join(REPO, "build", "mac", "release", "LKina")
BENCH = os.path.join(REPO, "benchmarks")
PREP = Path("/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks/redock_benchmark/prepared")
META = Path("/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks/redock_benchmark/results")

# system -> (metal_mode, tag) ; receptor/ligand from the redock prepared set
SYSTEMS = [
    ("4JC_ZN", "zn", "real_systems/4JC_core_rec.pdbqt", "4JC305_lig.pdbqt", (-2.55, 2.29, 85.60)),
    ("1A8E_FE", "fe3", "redock_benchmark/prepared/1A8E_FE_rec.pdbqt", "redock_benchmark/prepared/1A8E_FE_lig.pdbqt", None),
    ("1A2V_CU", "cu2_jt", "redock_benchmark/prepared/1A2V_CU_rec.pdbqt", "redock_benchmark/prepared/1A2V_CU_lig.pdbqt", None),
]
WEIGHTS = ["0.0", "0.1", "0.3", "0.5"]


def run(cmd, timeout=900):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"


def lig_center_from_meta(tag):
    f = META / f"{tag}_meta.json"
    if f.exists():
        d = json.load(open(f))
        return tuple(d["lig_center"]), [max(s + 6.0, 14) for s in d["span"]]
    return None


def donor_metal_dist(docked_pdbqt, rec_pdbqt):
    """Min distance from any ligand N/O/S heavy atom to any receptor metal.
    Atom types are case-insensitive: obabel writes Zn/Fe/Cu (mixed case),
    LKina AD4 writes ZN/FE/CU; ligand donors N/NA/NS/O/OA/OS/S/SA."""
    metals = []
    for ln in open(rec_pdbqt, errors="ignore"):
        if ln.startswith(("ATOM", "HETATM")):
            t = ln[77:79].strip()
            if t.lower() in ("zn", "fe", "cu", "mg", "ca", "mn", "co", "ni"):
                metals.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    if not metals:
        return None
    txt = open(docked_pdbqt, errors="ignore").read()
    seg = txt.split("MODEL")[1] if "MODEL" in txt else txt
    import math
    best = None
    for ln in seg.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            t = ln[77:79].strip()
            if t.lower() not in ("n", "na", "ns", "o", "oa", "os", "s", "sa"):
                continue
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            for m in metals:
                d = math.dist(xyz, m)
                best = d if best is None else min(best, d)
    return round(best, 2) if best is not None else None


def main():
    results = {}
    for name, mode, rec_rel, lig_rel, fixed_center in SYSTEMS:
        rec = PREP if rec_rel.startswith("redock") else Path("/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks")
        rec_p = Path("/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks") / rec_rel
        lig_p = Path("/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks") / lig_rel
        if not rec_p.exists() or not lig_p.exists():
            print(f"[{name}] missing {rec_p} or {lig_p}, skipping")
            continue
        tag = name.split("_")[0] + "_" + name.split("_")[1]
        if fixed_center:
            center, box = fixed_center, [25, 25, 25]
        else:
            mm = lig_center_from_meta(name)
            if not mm:
                print(f"[{name}] no meta, skipping"); continue
            center, box = mm
        print(f"== {name} (mode={mode}) center={center} box={box}")
        row = {"mode": mode, "w": {}}
        for w in WEIGHTS:
            out = f"/tmp/softw_{name.replace('_','')}_w{w.replace('.','')}.pdbqt"
            cmd = [LKINA, "--scoring", "ad4", "--generate_maps",
                   "--metal_mode", mode, "--metal_soft_weight", w,
                   "--receptor", str(rec_p), "--ligand", str(lig_p),
                   "--center_x", str(center[0]), "--center_y", str(center[1]),
                   "--center_z", str(center[2]),
                   "--size_x", str(box[0]), "--size_y", str(box[1]), "--size_z", str(box[2]),
                   "--exhaustiveness", "16", "--seed", "42", "--cpu", "6",
                   "--out", out, "--verbosity", "1"]
            rc, so, se = run(cmd)
            entry = {"rc": rc}
            if rc == 0 and os.path.exists(out):
                entry["donor_metal"] = donor_metal_dist(out, str(rec_p))
                m = re.search(r"^\s*1\s+(-?[0-9.]+)\s", so, re.M)
                entry["energy"] = float(m.group(1)) if m else None
                g = re.search(r"METAL_GEO_E:\s*(-?[0-9.]+)", open(out).read())
                entry["metal_geo_e"] = float(g.group(1)) if g else None
            else:
                entry["fail"] = (se or "").strip().splitlines()[-1][:100] if (se or "").strip() else "no output"
            row["w"][w] = entry
            print(f"  w={w:4s} rc={rc} donor-metal={entry.get('donor_metal')} "
                  f"E={entry.get('energy')} geo={entry.get('metal_geo_e')}")
        results[name] = row

    outpath = os.path.join(BENCH, "soft_weight_results.json")
    json.dump(results, open(outpath, "w"), indent=1)
    print("saved:", outpath)


if __name__ == "__main__":
    main()
