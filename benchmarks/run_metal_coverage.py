#!/usr/bin/env python3
"""Full metal-coverage benchmark for LKina.

Protocol for each metal_mode:
  A) score_only at ideal coordination geometry  -> METAL_COORD distance check
     (proves the engine's coordination model recognizes the geometry)
  B) global docking from displaced probe (5 A) -> recovery of coordination
     (proves the search actually recovers M-L coordination)
  C) Vina 1.2.7 baseline on the SAME synthetic system -> no coordination model
"""
import os, sys, math, subprocess, re, json

BENCH = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks"
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
VINA  = "/Users/luoxiaowen/Desktop/LKDock/源文件/LKDock_v4.0_Mac/engine/vina127"
SYN   = os.path.join(BENCH, "synthetic_metals")
os.makedirs(SYN, exist_ok=True)

MH_DIRS = [(1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0),(0,0,-1)]
TZ_DIRS = [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
SQ_DIRS = [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
LIN_DIRS= [(1,0,0),(-1,0,0)]

GEOM = {
  "zn":("TZ",2.00),"mg":("MH",2.05),"ca":("MH",2.35),"mn":("MH",2.20),
  "fe":("MH",2.10),"co":("MH",2.05),"ni":("SQ",1.95),"cu":("SQ",2.00),
  "pt":("SQ",2.05),"pd":("SQ",2.05),"ru":("MH",2.10),"ir":("MH",2.10),
  "au":("LIN",2.30),"rh":("MH",2.10),"ag":("LIN",2.30),"tc":("MH",2.10),
  "re":("MH",2.10),"os":("MH",2.10),"cd":("TZ",2.30),"hg":("LIN",2.35),
  "tl":("MH",2.40),"pb":("MH",2.45),"sb":("MH",2.30),"bi":("MH",2.40),
  "as":("MH",2.30),"na":("MH",2.40),"k":("MH",2.70),"li":("MH",2.10),
  "al":("MH",2.00),"sr":("MH",2.50),"ba":("MH",2.70),"v":("MH",2.00),
  "cr":("MH",2.00),"ti":("MH",2.05),"sc":("MH",2.10),"y":("MH",2.25),
  "zr":("MH",2.15),"nb":("MH",2.10),"hf":("MH",2.15),"ta":("MH",2.10),
  "w":("MH",2.05),"mo":("MH",2.05),"la":("MH",2.50),"ce":("MH",2.45),
  "pr":("MH",2.40),"nd":("MH",2.40),"sm":("MH",2.35),"eu":("MH",2.40),
  "gd":("MH",2.35),"tb":("MH",2.35),"dy":("MH",2.35),"ho":("MH",2.30),
  "er":("MH",2.30),"tm":("MH",2.30),"yb":("MH",2.30),"lu":("MH",2.25),
  "se":("MH",2.30),"ge":("MH",2.10),"ga":("MH",2.05),"in":("MH",2.25),
  "sn":("MH",2.30),
}
GEOM["cu2_jt"]=("JT",2.00); GEOM["mn3_jt"]=("JT",2.10)

TOKEN = {
  "zn":"Zn","mg":"Mg","ca":"Ca","mn":"Mn","fe":"Fe","co":"Co","ni":"Ni","cu":"Cu",
  "pt":"Pt","pd":"Pd","ru":"Ru","ir":"Ir","au":"Au","rh":"Rh","ag":"Ag","tc":"Tc",
  "re":"Re","os":"Os","cd":"Cd","hg":"Hg","tl":"Tl","pb":"Pb","sb":"Sb","bi":"Bi",
  "as":"As","na":"Na","k":"K","li":"Li","al":"Al","sr":"Sr","ba":"Ba","v":"V",
  "cr":"Cr","ti":"Ti","sc":"Sc","y":"Y","zr":"Zr","nb":"Nb","hf":"Hf","ta":"Ta",
  "w":"W","mo":"Mo","la":"La","ce":"Ce","pr":"Pr","nd":"Nd","sm":"Sm","eu":"Eu",
  "gd":"Gd","tb":"Tb","dy":"Dy","ho":"Ho","er":"Er","tm":"Tm","yb":"Yb","lu":"Lu",
  "se":"Se","ge":"Ge","ga":"Ga","in":"In","sn":"Sn","cu2_jt":"Cu","mn3_jt":"Mn",
}

def norm(v):
    l = math.sqrt(sum(c*c for c in v))
    return tuple(c/l for c in v)

def _atom_line(record, serial, name, resn, chain, resi, x, y, z, q, typ):
    s = f"{record:<6s}{serial:5d} "
    s += f"{name:>4s} "
    s += f"{resn:>3s}"
    s += f" {chain:1s}{resi:4d}"
    s += f"    {x:8.3f}{y:8.3f}{z:8.3f}"
    s += f"{1.00:6.2f}{0.00:6.2f}"
    s += f"    {q:>6s}"
    s += f"  {typ:<2s}"
    return s

def write_receptor(path, metal_token, donors):
    lines = []
    serial = 1
    for i, d in enumerate(donors):
        typ = "NA" if i % 2 == 0 else "OA"
        name = "ND1" if typ == "NA" else "OD1"
        resn = "HIS" if typ == "NA" else "ASP"
        lines.append(_atom_line("HETATM", serial, name, resn, "A", 301, *d, "+0.000", typ))
        serial += 1
    lines.append(_atom_line("HETATM", serial, "M", metal_token, "A", 302, 0.0, 0.0, 0.0, "+0.000", metal_token))
    lines.append("END")
    open(path, "w").write("\n".join(lines) + "\n")

def write_ligand(path, probe_z, probe_x=0.0):
    lines = [
      "REMARK  probe ligand (NA donor)",
      "ROOT",
      _atom_line("ATOM", 1, "N", "LIG", "A", 1, probe_x, 0.0, probe_z, "-0.350", "NA"),
      _atom_line("ATOM", 2, "C", "LIG", "A", 1, probe_x, 0.0, probe_z - 1.45, "+0.100", "C"),
      _atom_line("ATOM", 3, "H", "LIG", "A", 1, probe_x + 0.9, 0.0, probe_z - 1.45, "+0.100", "HD"),
      "ENDROOT",
      "TORSDOF 0",
    ]
    open(path, "w").write("\n".join(lines) + "\n")

def run(cmd, timeout=180):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"

def metal_coord_from_output(stdout):
    """Extract METAL_COORD distance from LKina verbose/score output."""
    m = re.findall(r"METAL_COORD[:\s]+(\w+)\s+donor=(\w+)\s+d=([0-9.]+)\s*A", stdout)
    return m

def donor_dist_from_outfile(path, metal_xyz=(0,0,0)):
    """Nearest ligand-atom distance to metal (origin) in the 1st model."""
    best = 99.0
    in_model = False
    for l in open(path):
        if l.startswith("MODEL"): in_model = True
        elif l.startswith("ENDMDL"): break
        elif in_model and l.startswith(("ATOM","HETATM")):
            try:
                d = math.sqrt(float(l[30:38])**2 + float(l[38:46])**2 + float(l[46:54])**2)
                if d < best: best = d
            except: pass
    return best

def main():
    results = []
    for mode in sorted(GEOM.keys()):
        cls, d0 = GEOM[mode]
        tok = TOKEN[mode]
        if cls == "TZ":   donors = [norm(v) for v in TZ_DIRS[:2]]
        elif cls == "SQ": donors = [norm(v) for v in SQ_DIRS[:2]]
        elif cls == "LIN":donors = [norm(v) for v in LIN_DIRS]
        else:             donors = [norm(v) for v in MH_DIRS[:2]]
        donors = [tuple(c*d0 for c in v) for v in donors]

        d = os.path.join(SYN, mode); os.makedirs(d, exist_ok=True)
        rec = os.path.join(d, "rec.pdbqt")
        lig_ideal = os.path.join(d, "lig_ideal.pdbqt")   # N at +z d0 (ideal)
        lig_disp  = os.path.join(d, "lig_disp.pdbqt")    # N at +z d0+5 (displaced)
        write_receptor(rec, tok, donors)
        write_ligand(lig_ideal, d0)
        write_ligand(lig_disp, d0 + 5.0)
        box = (0, 0, 0.5, 14, 14, 14)

        # --- A: score_only at ideal geometry ---
        rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", mode,
                          "--receptor", rec, "--ligand", lig_ideal, "--score_only",
                          "--center_x", "0", "--center_y", "0", "--center_z", "0.5",
                          "--size_x", "14", "--size_y", "14", "--size_z", "14",
                          "--verbosity", "1", "--cpu", "4"])
        coord_ok = False; coord_d = None
        if rc == 0:
            m = metal_coord_from_output(so + se)
            if m:
                coord_d = float(m[0][2]); coord_ok = abs(coord_d - d0) < 0.6

        # --- B: global docking from displaced probe ---
        out = os.path.join(d, "out.pdbqt")
        rc2, so2, se2 = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", mode,
                             "--receptor", rec, "--ligand", lig_disp,
                             "--center_x", "0", "--center_y", "0", "--center_z", "0.5",
                             "--size_x", "14", "--size_y", "14", "--size_z", "14",
                             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                             "--out", out, "--verbosity", "0"])
        dock_ok = rc2 == 0 and os.path.exists(out) and os.path.getsize(out) > 100
        d_lig = donor_dist_from_outfile(out) if dock_ok else None
        recover = dock_ok and d_lig is not None and d_lig < d0 + 0.8

        # --- C: Vina 1.2.7 baseline (same system, vina scoring) ---
        vout = os.path.join(d, "vina_out.pdbqt")
        rc3, so3, se3 = run([VINA, "--receptor", rec, "--ligand", lig_disp,
                             "--center_x", "0", "--center_y", "0", "--center_z", "0.5",
                             "--size_x", "14", "--size_y", "14", "--size_z", "14",
                             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                             "--out", vout, "--verbosity", "0"])
        vina_ok = rc3 == 0 and os.path.exists(vout) and os.path.getsize(vout) > 100
        vina_d = donor_dist_from_outfile(vout) if vina_ok else None
        vina_coord = vina_ok and vina_d is not None and vina_d < d0 + 0.8

        results.append({
            "mode": mode, "token": tok, "pseudo": cls, "d0": d0,
            "score_ok": rc == 0, "coord_ok": coord_ok, "coord_d": coord_d,
            "dock_ok": dock_ok, "d_lig": d_lig, "recover": recover,
            "vina_ok": vina_ok, "vina_d": vina_d, "vina_coord": vina_coord,
        })
        print(f"{mode:7s} {tok:3s} {cls:3s} d0={d0:.2f} | A:coord={coord_d if coord_d else 'NA':<5} "
              f"({coord_ok and 'OK' or 'x'}) | B:dock={d_lig if d_lig else 99:.2f} "
              f"({recover and 'RECOVER' or 'x'}) | C:vina={vina_d if vina_d else 99:.2f} "
              f"({vina_coord and 'coord' or 'no-coord'})")

    with open(os.path.join(BENCH, "metal_coverage_results.json"), "w") as f:
        json.dump(results, f, indent=2)

    n = len(results)
    n_coord = sum(1 for r in results if r["coord_ok"])
    n_recover = sum(1 for r in results if r["recover"])
    n_vina = sum(1 for r in results if r["vina_coord"])
    print(f"\n=== SUMMARY: {n} metal modes | geometry-recognition {n_coord}/{n} | "
          f"docking-recovery {n_recover}/{n} | vina-coordination {n_vina}/{n} ===")

if __name__ == "__main__":
    main()
