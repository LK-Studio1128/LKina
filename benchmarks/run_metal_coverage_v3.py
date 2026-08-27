#!/usr/bin/env python3
"""FINAL metal-coverage benchmark (v3).

Per metal_mode (63 modes):
  - global docking with LKina (metal_mode explicit) from displaced probe
  - record: success, METAL_COORD donor/distance (engine's own report),
            nearest donor-N distance to metal, energy
  - Vina 1.2.7 baseline on identical synthetic system
  - score_only energy at ideal vs far (coordination potential check)
"""
import os, sys, math, subprocess, re, json

BENCH = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks"
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
VINA  = "/Users/luoxiaowen/Desktop/LKDock/源文件/LKDock_v4.0_Mac/engine/vina127"
SYN   = os.path.join(BENCH, "synthetic_metals_v3")
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

def build(mode, d0, full):
    enter = norm(full[-1])
    rec_donors = [tuple(c*d0 for c in norm(v)) for v in full[:-1]]
    return rec_donors, enter

def write_receptor(path, tok, donors):
    lines = []
    s = 1
    lines.append(_atom_line("ATOM", s, "CA", "ALA", "A", 301, -3.0, 0.0, 0.0, "+0.100", "C")); s += 1
    for i, dd in enumerate(donors):
        typ = "NA" if i % 2 == 0 else "OA"
        name = "ND1" if typ == "NA" else "OD1"
        resn = "HIS" if typ == "NA" else "ASP"
        lines.append(_atom_line("ATOM", s, name, resn, "A", 301, *dd, "+0.000", typ)); s += 1
    lines.append(_atom_line("HETATM", s, "M", tok, "A", 302, 0.0, 0.0, 0.0, "+0.000", tok))
    lines.append("END")
    open(path, "w").write("\n".join(lines) + "\n")

def write_ligand(path, enter, dist):
    nx, ny, nz = enter[0]*dist, enter[1]*dist, enter[2]*dist
    open(path, "w").write("\n".join([
      "REMARK  probe ligand (NA donor)",
      "ROOT",
      _atom_line("ATOM", 1, "N", "LIG", "A", 1, nx, ny, nz, "-0.350", "NA"),
      _atom_line("ATOM", 2, "C", "LIG", "A", 1, nx, ny, nz - 1.45, "+0.100", "C"),
      _atom_line("ATOM", 3, "H", "LIG", "A", 1, nx + 0.9, ny, nz - 1.45, "+0.100", "HD"),
      "ENDROOT", "TORSDOF 0"]) + "\n")

def run(cmd, timeout=180):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"

def first_model(path):
    """Return (list_of_atoms, remarks_dict) for the first MODEL (or all atoms
    if no MODEL tags, i.e. local_only style output)."""
    atoms, remarks, in_model, started = [], {}, False, False
    for l in open(path):
        if l.startswith("MODEL"):
            atoms, remarks, in_model, started = [], {}, True, True
        elif l.startswith("ENDMDL"):
            if started: break
        elif l.startswith("REMARK"):
            m = re.match(r"REMARK\s+(\w[\w_]*):\s+(.+)", l)
            if m: remarks[m.group(1)] = m.group(2)
        elif l.startswith(("ATOM", "HETATM")):
            try:
                atoms.append((l[12:16].strip(), float(l[30:38]), float(l[38:46]),
                              float(l[46:54]), l[77:79].strip()))
            except: pass
    return atoms, remarks

def dist_to_metal(atoms):
    best = 99.0
    for a in atoms:
        d = math.sqrt(a[1]**2 + a[2]**2 + a[3]**2)
        if d < best: best = d
    return best

def main():
    results = []
    for mode in sorted(GEOM.keys()):
        cls, d0 = GEOM[mode]
        tok = TOKEN[mode]
        if cls == "TZ":   full = TZ_DIRS
        elif cls == "SQ": full = SQ_DIRS
        elif cls == "LIN":full = LIN_DIRS
        else:             full = MH_DIRS

        rec_donors, enter = build(mode, d0, full)
        d = os.path.join(SYN, mode); os.makedirs(d, exist_ok=True)
        rec = os.path.join(d, "rec.pdbqt")
        lig_ideal = os.path.join(d, "lig_ideal.pdbqt")
        lig_disp  = os.path.join(d, "lig_disp.pdbqt")
        write_receptor(rec, tok, rec_donors)
        write_ligand(lig_ideal, enter, d0)
        write_ligand(lig_disp, enter, d0 + 5.0)
        box = (0, 0, 0, 12, 12, 12)

        # energy check: ideal vs far (coordination potential)
        rc0, so0, _ = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", mode,
                           "--receptor", rec, "--ligand", lig_ideal, "--score_only",
                           "--center_x", "0", "--center_y", "0", "--center_z", "0",
                           "--size_x", "12", "--size_y", "12", "--size_z", "12",
                           "--verbosity", "0", "--cpu", "4"])
        e_ideal = None
        if rc0 == 0:
            m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[0-9.]+)", so0)
            if m: e_ideal = float(m.group(1))
        rc0b, so0b, _ = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", mode,
                             "--receptor", rec, "--ligand", lig_disp, "--score_only",
                             "--center_x", "0", "--center_y", "0", "--center_z", "0",
                             "--size_x", "12", "--size_y", "12", "--size_z", "12",
                             "--verbosity", "0", "--cpu", "4"])
        e_far = None
        if rc0b == 0:
            m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[0-9.]+)", so0b)
            if m: e_far = float(m.group(1))

        # LKina global docking
        out = os.path.join(d, "lkina_out.pdbqt")
        rc1, so1, se1 = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", mode,
                             "--receptor", rec, "--ligand", lig_disp,
                             "--center_x", "0", "--center_y", "0", "--center_z", "0",
                             "--size_x", "12", "--size_y", "12", "--size_z", "12",
                             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                             "--out", out, "--verbosity", "0"])
        atoms, remarks = first_model(out) if (rc1 == 0 and os.path.exists(out) and os.path.getsize(out) > 100) else ([], {})
        d_lig = dist_to_metal(atoms) if atoms else None
        coord = remarks.get("METAL_COORD", "")

        # Vina baseline
        vout = os.path.join(d, "vina_out.pdbqt")
        rc2, so2, se2 = run([VINA, "--receptor", rec, "--ligand", lig_disp,
                             "--center_x", "0", "--center_y", "0", "--center_z", "0",
                             "--size_x", "12", "--size_y", "12", "--size_z", "12",
                             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                             "--out", vout, "--verbosity", "0"])
        vatoms, _ = first_model(vout) if (rc2 == 0 and os.path.exists(vout) and os.path.getsize(vout) > 100) else ([], {})
        vina_d = dist_to_metal(vatoms) if vatoms else None

        results.append({
            "mode": mode, "token": tok, "pseudo": cls, "d0": d0,
            "e_ideal": e_ideal, "e_far": e_far,
            "lkina_rc": rc1, "d_lig": d_lig, "coord": coord,
            "vina_rc": rc2, "vina_d": vina_d,
            "err": (se1 or so1)[-120:] if rc1 != 0 else "",
            "verr": (se2 or so2)[-120:] if rc2 != 0 else "",
        })
        print(f"{mode:7s} {tok:3s} {cls:3s} d0={d0:.2f} | E_ideal={e_ideal if e_ideal is not None else 99:6.1f} "
              f"E_far={e_far if e_far is not None else 99:6.1f} | LKina d={d_lig if d_lig else 99:4.2f} "
              f"{'OK' if rc1==0 else 'FAIL'} | Vina d={vina_d if vina_d else 99:4.2f} "
              f"{'OK' if rc2==0 else 'FAIL'} | {coord[:32]}")

    with open(os.path.join(BENCH, "metal_coverage_results_v3.json"), "w") as f:
        json.dump(results, f, indent=2)

    n = len(results)
    n_ok = sum(1 for r in results if r["lkina_rc"] == 0)
    n_vok = sum(1 for r in results if r["vina_rc"] == 0)
    # coordination potential: e_ideal < e_far means the metal well exists
    n_well = sum(1 for r in results if r["e_ideal"] is not None and r["e_far"] is not None
                 and r["e_ideal"] < r["e_far"] + 0.3)
    print(f"\n=== SUMMARY: {n} modes | LKina-dock-OK {n_ok}/{n} | Vina-OK {n_vok}/{n} | "
          f"metal-well(E_ideal<E_far) {n_well}/{n} ===")

if __name__ == "__main__":
    main()
