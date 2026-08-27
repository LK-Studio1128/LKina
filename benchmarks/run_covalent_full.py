#!/usr/bin/env python3
"""Full covalent framework test: P1-P4 tiers, NAC detection, C3 two-step,
6 presets, reactive energy-well scan, and inline generate_maps probe coverage.
"""
import os, sys, math, subprocess, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pdbqt_util import atom_line, norm

BENCH = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks"
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
RXN   = os.path.join(BENCH, "covalent_full")
os.makedirs(RXN, exist_ok=True)

PRESETS = {
  "cys_michael":   {"rec": ("SG", "CYS", "S"), "lig": "C",  "r0": 1.82, "theta": 109.5},
  "cys_sn2":       {"rec": ("SG", "CYS", "S"), "lig": "C",  "r0": 1.82, "theta": 180.0},
  "ser_covalent":  {"rec": ("OG", "SER", "O"), "lig": "C",  "r0": 1.34, "theta": 109.5},
  "lys_targeting": {"rec": ("NZ", "LYS", "N"), "lig": "C",  "r0": 1.47, "theta": 109.5},
  "boronic_acid":  {"rec": ("OG", "SER", "O"), "lig": "B",  "r0": 1.47, "theta": None},
  "tyr_covalent":  {"rec": ("OH", "TYR", "O"), "lig": "C",  "r0": 1.38, "theta": 109.5},
}

def build_system(name, p):
    d = os.path.join(RXN, name); os.makedirs(d, exist_ok=True)
    rec_atom_name, resn, rec_typ = p["rec"]
    r0 = p["r0"]
    rec_lines = []
    rec_lines.append(atom_line("ATOM", 1, "CA", "ALA", "A", 1, -2.0, 0.0, 0.0, "+0.100", "C"))
    rec_lines.append(atom_line("ATOM", 2, "CB", "ALA", "A", 1, -1.0, 0.0, 0.0, "+0.100", "C"))
    rec_lines.append(atom_line("ATOM", 3, rec_atom_name, resn, "A", 1, 0.0, 0.0, 0.0, "-0.200", rec_typ))
    rec_lines.append(atom_line("ATOM", 4, "CB", resn, "A", 1, 1.8, 0.0, 0.0, "+0.100", "C"))
    rec_lines.append("END")
    rec_path = os.path.join(d, "rec.pdbqt")
    open(rec_path, "w").write("\n".join(rec_lines) + "\n")
    lig_lines = [
      "REMARK  reactive probe ligand",
      "ROOT",
      atom_line("ATOM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, r0, "-0.100", p["lig"]),
      atom_line("ATOM", 2, "C2", "LIG", "A", 1, 0.0, 0.0, -1.5, "+0.100", "C"),
      atom_line("ATOM", 3, "H1", "LIG", "A", 1, 0.9, 0.0, r0, "+0.050", "HD"),
      "ENDROOT",
      "TORSDOF 0",
    ]
    lig_path = os.path.join(d, "lig.pdbqt")
    open(lig_path, "w").write("\n".join(lig_lines) + "\n")
    return d, rec_path, lig_path

def run(cmd, timeout=180):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"

def collect(path, extra_txt=""):
    txt = extra_txt
    if os.path.exists(path):
        try: txt += open(path).read()
        except: pass
    out = {}
    for l in txt.splitlines():
        for key in ("REACTIVE_NAC", "REACTIVE_DIST", "REACTIVE_ANGLE", "REACTIVE_DIST_E",
                    "REACTIVE_ANGLE_E", "VINA RESULT"):
            m = re.search(rf"{key}[:=]\s*([YESNO0-9.eE+-]+)", l)
            if m:
                v = m.group(1)
                if v in ("YES", "NO"): out[key] = v
                else:
                    try: out[key] = float(v)
                    except: pass
                break
    return out

BASE_ARGS = ["--scoring", "ad4", "--generate_maps",
             "--center_x", "0", "--center_y", "0", "--center_z", "0",
             "--size_x", "16", "--size_y", "16", "--size_z", "16",
             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4", "--verbosity", "0"]

results = {}
preset_rows = []

for name, p in PRESETS.items():
    d, rec, lig = build_system(name, p)
    anchor = "A:1:" + p["rec"][0]
    frame = "A:1:CB"
    rows = {}

    # ---- P1: distance-only ----
    out1 = os.path.join(d, "p1.pdbqt")
    rc1, so1, se1 = run([LKINA] + BASE_ARGS + ["--receptor", rec, "--ligand", lig,
        "--reactive_preset", name, "--reactive_rec_atom", anchor,
        "--reactive_lig_atom", "index:1", "--reactive_mode", "distance",
        "--out", out1])
    r1 = collect(out1, so1 + se1) if rc1 == 0 else {}
    rows["p1"] = {"rc": rc1, "nac": r1.get("REACTIVE_NAC"), "dist": r1.get("REACTIVE_DIST")}

    # ---- P1+P2: distance + angle (hybrid) ----
    out2 = os.path.join(d, "p12.pdbqt")
    rc2, so2, se2 = run([LKINA] + BASE_ARGS + ["--receptor", rec, "--ligand", lig,
        "--reactive_preset", name, "--reactive_rec_atom", anchor,
        "--reactive_lig_atom", "index:1", "--reactive_mode", "hybrid",
        "--reactive_frame_atom", frame, "--reactive_lig_frame_atom", "index:2",
        "--out", out2])
    r2 = collect(out2, so2 + se2) if rc2 == 0 else {}
    rows["p12"] = {"rc": rc2, "nac": r2.get("REACTIVE_NAC"), "dist": r2.get("REACTIVE_DIST"),
                   "angle": r2.get("REACTIVE_ANGLE")}

    # ---- P4: hybrid vdW scale sweep (0 / 0.5 / 1.0) ----
    p4 = {}
    for scale in ("0.0", "0.5", "1.0"):
        out4 = os.path.join(d, f"p4_vdw{scale.replace('.','')}.pdbqt")
        rc4, so4, se4 = run([LKINA] + BASE_ARGS + ["--receptor", rec, "--ligand", lig,
            "--reactive_preset", name, "--reactive_rec_atom", anchor,
            "--reactive_lig_atom", "index:1", "--reactive_mode", "hybrid",
            "--reactive_frame_atom", frame, "--reactive_lig_frame_atom", "index:2",
            "--reactive_hybrid_vdw_scale", scale,
            "--out", out4])
        r4 = collect(out4, so4 + se4) if rc4 == 0 else {}
        e = r4.get("VINA RESULT")
        p4[scale] = {"rc": rc4, "energy": e}
    rows["p4_vdw_scale"] = p4

    # ---- C3 two-step ----
    out3 = os.path.join(d, "two_step.pdbqt")
    rc3, so3, se3 = run([LKINA] + BASE_ARGS + ["--receptor", rec, "--ligand", lig,
        "--reactive_preset", name, "--reactive_rec_atom", anchor,
        "--reactive_lig_atom", "index:1", "--reactive_mode", "hybrid",
        "--reactive_frame_atom", frame, "--reactive_lig_frame_atom", "index:2",
        "--reactive_two_step", "--reactive_presample_dist", "10.0",
        "--out", out3])
    rows["c3_two_step"] = {"rc": rc3, "size": os.path.getsize(out3) if rc3 == 0 and os.path.exists(out3) else 0}

    preset_rows.append({"preset": name, "r0": p["r0"], "theta": p["theta"], **rows})
    print(f"{name:14s} P1 rc={rc1} nac={rows['p1']['nac']} | P12 rc={rc2} nac={rows['p12']['nac']} "
          f"angle={rows['p12']['angle']} | vdw0={p4['0.0']['energy']} vdw1={p4['1.0']['energy']} "
          f"| 2step rc={rc3} size={rows['c3_two_step']['size']}")

results["presets"] = preset_rows

# ---- Reactive energy-well scan (cys_michael: distance vs E) ----
scan = {"preset": "cys_michael", "points": []}
d0 = PRESETS["cys_michael"]["r0"]
d = os.path.join(RXN, "cys_michael")
rec = os.path.join(d, "rec.pdbqt")
for delta in [-0.8, -0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5, 0.8]:
    dist = d0 + delta
    lig = os.path.join(d, f"lig_scan_{dist:.1f}.pdbqt")
    open(lig, "w").write("\n".join([
      "REMARK", "ROOT",
      atom_line("ATOM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, dist, "-0.100", "C"),
      atom_line("ATOM", 2, "C2", "LIG", "A", 1, 0.0, 0.0, -1.5, "+0.100", "C"),
      atom_line("ATOM", 3, "H1", "LIG", "A", 1, 0.9, 0.0, dist, "+0.050", "HD"),
      "ENDROOT", "TORSDOF 0"]) + "\n")
    rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps",
        "--receptor", rec, "--ligand", lig, "--score_only",
        "--center_x", "0", "--center_y", "0", "--center_z", "0",
        "--size_x", "16", "--size_y", "16", "--size_z", "16",
        "--reactive_preset", "cys_michael", "--reactive_rec_atom", "A:1:SG",
        "--reactive_lig_atom", "index:1", "--reactive_mode", "distance",
        "--reactive_debug_energy", "--verbosity", "1", "--cpu", "4"])
    e_total = None; e_react = None
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[0-9.]+)", so)
    if m: e_total = float(m.group(1))
    m = re.search(r"Reactive Distance Energy\s*:\s*(-?[0-9.]+)", so + se)
    if m: e_react = float(m.group(1))
    scan["points"].append({"delta": delta, "dist": dist, "e_total": e_total, "e_reactive": e_react})
    print(f"  scan d={dist:.1f} E_total={e_total if e_total is not None else 99:7.2f} E_reactive={e_react if e_react is not None else 99:7.2f}")
results["well_scan"] = scan

# ---- generate_maps probe coverage (complex organic ligand) ----
gm = os.path.join(RXN, "genmaps"); os.makedirs(gm, exist_ok=True)
# a diverse organic ligand: C, N, OA, SA, HD, F, Cl, P
lig = os.path.join(gm, "lig.pdbqt")
open(lig, "w").write("\n".join([
  "REMARK  diverse organic ligand", "ROOT",
  atom_line("ATOM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, 0.0, "+0.100", "C"),
  atom_line("ATOM", 2, "N2", "LIG", "A", 1, 1.4, 0.0, 0.0, "-0.350", "NA"),
  atom_line("ATOM", 3, "O3", "LIG", "A", 1, 0.0, 1.4, 0.0, "-0.350", "OA"),
  atom_line("ATOM", 4, "S4", "LIG", "A", 1, 0.0, 0.0, 1.4, "-0.200", "SA"),
  atom_line("ATOM", 5, "F5", "LIG", "A", 1, 2.0, 0.5, 0.5, "-0.150", "F"),
  atom_line("ATOM", 6, "CL6", "LIG", "A", 1, -1.4, 0.0, 0.0, "-0.150", "Cl"),
  atom_line("ATOM", 7, "P7", "LIG", "A", 1, 0.0, -1.4, 0.0, "+0.500", "P"),
  atom_line("ATOM", 8, "H8", "LIG", "A", 1, 2.2, 0.0, 0.0, "+0.100", "HD"),
  "ENDROOT", "TORSDOF 0"]) + "\n")
rec = os.path.join(gm, "rec.pdbqt")
open(rec, "w").write("\n".join([
  atom_line("ATOM", 1, "CA", "ALA", "A", 301, -3.0, 0.0, 0.0, "+0.100", "C"),
  atom_line("ATOM", 2, "ND1", "HIS", "A", 301, 3.0, 0.0, 0.0, "+0.000", "NA"),
  atom_line("ATOM", 3, "OD1", "ASP", "A", 301, 0.0, 3.0, 0.0, "+0.000", "OA"),
  atom_line("ATOM", 4, "SG", "CYS", "A", 301, 0.0, 0.0, 3.0, "+0.000", "SA"),
  "END"]) + "\n")
mapsdir = os.path.join(gm, "maps")
os.makedirs(mapsdir, exist_ok=True)
rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps",
    "--receptor", rec, "--ligand", lig,
    "--center_x", "0", "--center_y", "0", "--center_z", "0",
    "--size_x", "14", "--size_y", "14", "--size_z", "14",
    "--write_maps", os.path.join(mapsdir, "g"), "--verbosity", "1", "--cpu", "4"])
maps = sorted(f for f in os.listdir(mapsdir) if f.endswith(".map"))
results["generate_maps"] = {"rc": rc, "n_maps": len(maps), "maps": maps,
                            "probe_types": ["C","NA","OA","SA","F","Cl","P","HD"]}
print(f"\ngenerate_maps: rc={rc} maps={len(maps)} -> {maps}")

with open(os.path.join(BENCH, "covalent_full_results.json"), "w") as f:
    json.dump(results, f, indent=2)

n_p1 = sum(1 for r in preset_rows if r["p1"]["rc"] == 0)
n_p12 = sum(1 for r in preset_rows if r["p12"]["rc"] == 0)
n_p4 = sum(1 for r in preset_rows if all(r["p4_vdw_scale"][k]["rc"] == 0 for k in ("0.0","0.5","1.0")))
n_c3 = sum(1 for r in preset_rows if r["c3_two_step"]["rc"] == 0)
n_nac_yes = sum(1 for r in preset_rows if r["p12"].get("nac") == "YES")
print(f"\n=== SUMMARY: P1 {n_p1}/6 | P1+P2 {n_p12}/6 | P4 sweep {n_p4}/6 | C3 {n_c3}/6 | NAC=YES {n_nac_yes}/6 | "
      f"well-scan {len(scan['points'])} pts | maps {results['generate_maps']['n_maps']} ===")
