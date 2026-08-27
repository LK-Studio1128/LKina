#!/usr/bin/env python3
"""Feature-family tests: pseudoatom geometry, BVS oxidation-state inference,
Jahn-Teller, water bridge, and metal-as-ligand geometry QC.

Each family produces a JSON row so the paper can report per-feature evidence.
"""
import os, sys, math, subprocess, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pdbqt_util import atom_line, norm

BENCH = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks"
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
OUT   = os.path.join(BENCH, "feature_family_tests")
os.makedirs(OUT, exist_ok=True)

def write_simple_receptor(path, tok, donors):
    """Metal at origin, donors at given positions, 1 C backbone atom."""
    lines = []
    s = 1
    lines.append(atom_line("ATOM", s, "CA", "ALA", "A", 301, -3.0, 0.0, 0.0, "+0.100", "C")); s += 1
    for i, dd in enumerate(donors):
        typ = "NA" if i % 2 == 0 else "OA"
        name = "ND1" if typ == "NA" else "OD1"
        resn = "HIS" if typ == "NA" else "ASP"
        lines.append(atom_line("ATOM", s, name, resn, "A", 301, *dd, "+0.000", typ)); s += 1
    lines.append(atom_line("HETATM", s, "M", tok, "A", 302, 0.0, 0.0, 0.0, "+0.000", tok))
    lines.append("END")
    open(path, "w").write("\n".join(lines) + "\n")

def write_ligand(path, pos, typ="NA", name="N"):
    nx, ny, nz = pos
    open(path, "w").write("\n".join([
      "REMARK  probe",
      "ROOT",
      atom_line("ATOM", 1, name, "LIG", "A", 1, nx, ny, nz, "-0.350", typ),
      atom_line("ATOM", 2, "C", "LIG", "A", 1, nx, ny, nz - 1.45, "+0.100", "C"),
      atom_line("ATOM", 3, "H", "LIG", "A", 1, nx + 0.9, ny, nz - 1.45, "+0.100", "HD"),
      "ENDROOT", "TORSDOF 0"]) + "\n")

def run(cmd, timeout=180):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"

def parse_remarks(path):
    """Return dict of all REMARK key:value from first MODEL block."""
    rem = {}
    for l in open(path):
        if l.startswith("MODEL"): rem = {}
        elif l.startswith("ENDMDL"): break
        elif l.startswith("REMARK"):
            m = re.match(r"REMARK\s+(\w[\w_]*):\s*(.*)", l)
            if m: rem[m.group(1)] = m.group(2)
    return rem

results = {}

# ============================================================
# A. Pseudoatom geometry checks (TZ / SQ / MH / JT)
# ============================================================
geo_cases = [
    ("zn",      "TZ", 2.00, "Zn", [(1,1,1),(1,-1,-1),(-1,1,-1)]),  # 3 occ -> 1 TZ
    ("pt",      "SQ", 2.05, "Pt", [(1,0,0),(0,1,0),(-1,0,0)]),      # 3 occ -> 1 SQ
    ("fe",      "MH", 2.10, "Fe", [(1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0)]),  # 5 occ -> 1 MH
    ("cu2_jt",  "JT", 2.00, "Cu", [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]),         # JT equatorial
    ("mn3_jt",  "JT", 2.10, "Mn", [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]),
]
geo_rows = []
for mode, cls, d0, tok, dirs in geo_cases:
    donors = [tuple(c*d0 for c in norm(v)) for v in dirs]
    d = os.path.join(OUT, "geo_" + mode); os.makedirs(d, exist_ok=True)
    rec = os.path.join(d, "rec.pdbqt")
    lig = os.path.join(d, "lig.pdbqt")
    write_simple_receptor(rec, tok, donors)
    write_ligand(lig, (0, 0, d0 + 2.0))
    rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", mode,
                      "--metal_geometry_check", "--receptor", rec,
                      "--ligand", lig, "--score_only",
                      "--center_x", "0", "--center_y", "0", "--center_z", "0",
                      "--size_x", "12", "--size_y", "12", "--size_z", "12",
                      "--verbosity", "1", "--cpu", "4"])
    txt = so + se
    # parse geometry-check block (printf format: %-6s %-4s %-4s %-4s %8.3f %8.3f %+8.3f quality)
    dists, quals, cn, vacancies = [], [], None, None
    for l in txt.splitlines():
        m = re.match(r"\s*\S+\s+\S+\s+\S+\s+\S+\s+([0-9.]+)\s+([0-9.]+)\s+[+-][0-9.]+\s+(ideal|good|fair|poor)", l)
        if m:
            dists.append(float(m.group(1))); quals.append(m.group(3))
        m2 = re.search(r"CN_receptor=(\d+)\s*/\s*max_CN=(\d+)\s+vacancies=(\d+)", l)
        if m2: cn, mx, vacancies = int(m2.group(1)), int(m2.group(2)), int(m2.group(3))
    n_ideal = sum(1 for q in quals if q in ("ideal", "good"))
    geo_rows.append({"mode": mode, "pseudo": cls, "d0": d0,
                     "n_donors_checked": len(dists),
                     "n_ideal": n_ideal,
                     "dists": [round(x,3) for x in dists],
                     "vacancies": vacancies,
                     "rc": rc})
    print(f"[A] {mode:7s} {cls} d0={d0} donors={len(dists)} ideal={n_ideal}/{len(dists)} vac={vacancies} rc={rc}")
results["pseudoatom_geometry"] = geo_rows

# ============================================================
# B. BVS oxidation-state inference (auto-detection, no --metal_mode)
# ============================================================
bvs_cases = [
    # (tok, donors at dist, expected mode, label)
    ("Fe", 2.10, "fe2", "Fe-6N@2.10 -> Fe2+"),
    ("Fe", 1.95, "fe3", "Fe-6N@1.95 -> Fe3+"),
    ("Cu", 2.10, "cu1", "Cu-4N@2.10 -> Cu1+"),
    ("Cu", 1.95, "cu2", "Cu-4N@1.95 -> Cu2+"),
    ("Mn", 2.20, "mn2", "Mn-6N@2.20 -> Mn2+"),
    ("Mn", 2.00, "mn3", "Mn-6N@2.00 -> Mn3+"),
    ("Co", 2.15, "co2", "Co-6N@2.15 -> Co2+"),
    ("Co", 1.90, "co3", "Co-6N@1.90 -> Co3+"),
    ("V",  2.05, "v4",  "V-6N@2.05  -> V4+"),
    ("V",  1.88, "v5",  "V-6N@1.88  -> V5+"),
    ("Mo", 2.10, "mo4", "Mo-6O@2.10 -> Mo4+"),
    ("Mo", 1.90, "mo6", "Mo-6O@1.90 -> Mo6+"),
    ("Ni", 2.05, "ni2", "Ni-4N@2.05 -> Ni2+"),
    ("Ni", 1.73, "ni3", "Ni-4N@1.73 -> Ni3+"),
]
bvs_rows = []
for tok, d0, expect, label in bvs_cases:
    n = 6 if tok in ("Fe","Mn","Co","V","Mo") else 4
    dirs = []
    for k in range(n):
        # evenly spread donor directions (alternate axes for 4, octahedral for 6)
        if n == 6:
            v = [(1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0),(0,0,-1)][k]
        else:
            v = [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)][k]
        dirs.append(tuple(c*d0 for c in v))
    d = os.path.join(OUT, "bvs_" + label.replace(" ","_").replace("+","p").replace("@","_").replace("->","_to_").replace(":",""))
    os.makedirs(d, exist_ok=True)
    rec = os.path.join(d, "rec.pdbqt")
    # donors all OA for Mo (O-preference), NA for others
    lines = []
    s = 1
    lines.append(atom_line("ATOM", s, "CA", "ALA", "A", 301, -3.0, 0.0, 0.0, "+0.100", "C")); s += 1
    for dd in dirs:
        typ = "OA" if tok == "Mo" else "NA"
        name = "OD1" if typ == "OA" else "ND1"
        resn = "ASP" if typ == "OA" else "HIS"
        lines.append(atom_line("ATOM", s, name, resn, "A", 301, *dd, "+0.000", typ)); s += 1
    lines.append(atom_line("HETATM", s, "M", tok, "A", 302, 0.0, 0.0, 0.0, "+0.000", tok))
    lines.append("END")
    open(rec, "w").write("\n".join(lines) + "\n")
    lig = os.path.join(d, "lig.pdbqt")
    write_ligand(lig, (0, 0, d0 + 2.0))
    rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps",
                      "--receptor", rec, "--ligand", lig, "--score_only",
                      "--center_x", "0", "--center_y", "0", "--center_z", "0",
                      "--size_x", "12", "--size_y", "12", "--size_z", "12",
                      "--verbosity", "0", "--cpu", "4"])
    detected = None
    m = re.search(r"Auto-detected (\w+) metal in receptor", se or so)
    if m: detected = m.group(1)
    ok = detected == expect
    bvs_rows.append({"tok": tok, "d0": d0, "expected": expect, "detected": detected, "ok": ok, "label": label})
    print(f"[B] {tok:2s} @{d0:.2f} expected={expect:4s} detected={detected} {'OK' if ok else 'MISMATCH'}")
results["bvs_inference"] = bvs_rows

# ============================================================
# C. Water bridge (METAL_WATER_E) — Mg with 2 vacant octahedral sites
# ============================================================
wc = os.path.join(OUT, "waterbridge_mg"); os.makedirs(wc, exist_ok=True)
rec = os.path.join(wc, "rec.pdbqt")
# Mg octahedral max_CN=6, only 4 donors -> 2 vacant -> water sites
mh4 = [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
donors = [tuple(c*2.05 for c in v) for v in mh4]
lines = []
s = 1
lines.append(atom_line("ATOM", s, "CA", "ALA", "A", 301, -3.0, 0.0, 0.0, "+0.100", "C")); s += 1
for i, dd in enumerate(donors):
    typ = "NA" if i % 2 == 0 else "OA"
    name = "ND1" if typ == "NA" else "OD1"
    resn = "HIS" if typ == "NA" else "ASP"
    lines.append(atom_line("ATOM", s, name, resn, "A", 301, *dd, "+0.000", typ)); s += 1
lines.append(atom_line("HETATM", s, "M", "Mg", "A", 302, 0.0, 0.0, 0.0, "+0.000", "Mg"))
lines.append("END")
open(rec, "w").write("\n".join(lines) + "\n")
lig = os.path.join(wc, "lig.pdbqt")
write_ligand(lig, (0, 0, 2.05))
out = os.path.join(wc, "out.pdbqt")
# Global docking: full MODEL+REMARK output so METAL_RERANK / METAL_GEO_E /
# METAL_WATER_E / METAL_COORD are reported per pose.
rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", "mg",
                  "--receptor", rec, "--ligand", lig,
                  "--center_x", "0", "--center_y", "0", "--center_z", "0",
                  "--size_x", "12", "--size_y", "12", "--size_z", "12",
                  "--exhaustiveness", "16", "--seed", "42", "--cpu", "4",
                  "--out", out, "--verbosity", "0"])
rem = parse_remarks(out) if rc == 0 else {}
water_e = rem.get("METAL_WATER_E")
geo_e = rem.get("METAL_GEO_E")
coord = rem.get("METAL_COORD")
print(f"[C] water bridge: rc={rc} METAL_WATER_E={water_e} METAL_GEO_E={geo_e} COORD={coord}")
results["water_bridge"] = {"rc": rc, "metal_water_e": water_e, "metal_geo_e": geo_e,
                           "metal_coord": coord,
                           "n_vacant_sites": 2, "config": "Mg 4-donor (2 vacant octahedral)"}

# ============================================================
# D. Metal-as-ligand geometry QC — Pt(II) square-planar ligand
# ============================================================
def write_pt_ligand(path, mode):
    """Proper square-planar Pt(II) ligand in the xy plane.
    ideal:   4 ligands at 90 deg intervals (N-Pt-N = 90 deg)
    distorted: N-Pt-N = 60 deg (in-plane), Cl trans to the N bisector.
    """
    if mode == "ideal":
        n1, n2 = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
    else:  # distorted: N-Pt-N = 60 deg
        n1, n2 = (math.cos(math.radians(30)), math.sin(math.radians(30)), 0.0), \
                 (math.cos(math.radians(30)), -math.sin(math.radians(30)), 0.0)
    c1, c2 = (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0)
    d_pt_n = 2.05; d_pt_cl = 2.30
    lines = ["REMARK  Pt(II) sq-planar ligand", "ROOT"]
    lines.append(atom_line("ATOM", 1, "PT", "LIG", "A", 1, 0.0, 0.0, 0.0, "+0.500", "Pt"))
    for k, p in enumerate([n1, n2]):
        lines.append(atom_line("ATOM", 2+k, "N", "LIG", "A", 1, p[0]*d_pt_n, p[1]*d_pt_n, p[2]*d_pt_n, "-0.350", "NA"))
    for k, p in enumerate([c1, c2]):
        lines.append(atom_line("ATOM", 4+k, "CL", "LIG", "A", 1, p[0]*d_pt_cl, p[1]*d_pt_cl, p[2]*d_pt_cl, "-0.250", "Cl"))
    lines += ["ENDROOT", "TORSDOF 0"]
    open(path, "w").write("\n".join(lines) + "\n")

# receptor: a donor-rich "protein" surface above the ligand box (simple)
mc = os.path.join(OUT, "metal_as_ligand"); os.makedirs(mc, exist_ok=True)
rec = os.path.join(mc, "rec.pdbqt")
rec_lines = []
s = 1
for (rx, ry, rz, typ, name, resn) in [
    (2.0, 0.0, 2.0, "NA", "ND1", "HIS"), (0.0, 2.0, 2.0, "NA", "ND1", "HIS"),
    (-2.0, 0.0, 2.0, "OA", "OD1", "ASP"), (0.0, -2.0, 2.0, "OA", "OD1", "ASP"),
    (3.0, 0.0, -2.0, "OA", "OD1", "ASP"), (0.0, 3.0, -2.0, "NA", "ND1", "HIS"),
]:
    rec_lines.append(atom_line("ATOM", s, name, resn, "A", 301, rx, ry, rz, "+0.000", typ)); s += 1
rec_lines.append("END")
open(rec, "w").write("\n".join(rec_lines) + "\n")

lig_rows = []
for label, angle in [("ideal_sq_planar", 90.0), ("distorted", 60.0)]:
    lig = os.path.join(mc, f"lig_{label}.pdbqt")
    write_pt_ligand(lig, "ideal" if label.startswith("ideal") else "distorted")
    out = os.path.join(mc, f"out_{label}.pdbqt")
    rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps",
                      "--receptor", rec, "--ligand", lig,
                      "--center_x", "0", "--center_y", "0", "--center_z", "0",
                      "--size_x", "14", "--size_y", "14", "--size_z", "14",
                      "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                      "--ligand_metal_geometry_weight", "1.0",
                      "--out", out, "--verbosity", "0"])
    rem = parse_remarks(out) if rc == 0 else {}
    lm_geom = rem.get("LIGAND_METAL_GEOM")
    lm_rerank = rem.get("LIGAND_METAL_RERANK")
    lm_site = rem.get("LIGAND_METAL_SITE")
    print(f"[D] Pt ligand {label:16s} rc={rc} LIGAND_METAL_GEOM={lm_geom} RERANK={lm_rerank} SITE={lm_site}")
    lig_rows.append({"label": label, "angle": angle, "rc": rc,
                     "ligand_metal_geom": lm_geom, "ligand_metal_rerank": lm_rerank,
                     "ligand_metal_site": lm_site})
results["metal_as_ligand"] = lig_rows

with open(os.path.join(BENCH, "feature_family_results.json"), "w") as f:
    json.dump(results, f, indent=2)

# summary
nb = sum(1 for r in bvs_rows if r["ok"])
print(f"\n=== SUMMARY: geometry {sum(1 for r in geo_rows if r['rc']==0)}/{len(geo_rows)} | "
      f"BVS {nb}/{len(bvs_rows)} | water {results['water_bridge']['rc']==0} | "
      f"metal-ligand {sum(1 for r in lig_rows if r['rc']==0)}/{len(lig_rows)} ===")
