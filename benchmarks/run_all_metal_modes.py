#!/usr/bin/env python3
"""FULL metal-mode coverage benchmark (all 110 CLI metal modes).

For EVERY metal_mode registered by LKina v1.0.0:
  - synthetic coordination system (metal centre + donor shell + probe ligand)
  - LKina score_only at ideal vs far position (coordination-potential check)
  - LKina global docking  (success + METAL_COORD + nearest donor distance)
  - Vina 1.2.7 baseline    (shows which metals its XS type system can parse)
Output: metal_coverage_results_all.json  (per-mode rows)
"""
import os, sys, math, subprocess, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pdbqt_util import atom_line, norm, write_receptor, write_ligand

BENCH = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks"
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
VINA  = "/Users/luoxiaowen/Desktop/LKDock/源文件/LKDock_v4.0_Mac/engine/vina127"
SYN   = os.path.join(BENCH, "synthetic_metals_all_v3")
os.makedirs(SYN, exist_ok=True)

# pseudo class, ideal donor distance d0 (A), receptor atom token
# d0 taken from the mode's own r_eq overrides / Harding CSD survey values.
M = {}
def add(mode, cls, d0, tok): M[mode] = (cls, d0, tok)

# ---- M1 special: Zn tetrahedral ----
add("zn",      "TZ",  2.00, "Zn")
# ---- biological first-row TMs ----
add("mg",      "MH",  2.05, "Mg");  add("ca",   "MH",  2.35, "Ca")
add("mn",      "MH",  2.20, "Mn");  add("fe",   "MH",  2.10, "Fe")
add("co",      "MH",  2.05, "Co");  add("ni",   "SQ",  1.95, "Ni")
add("cu",      "SQ",  2.00, "Cu")
# ---- medicinal organometallics ----
add("pt",      "SQ",  2.05, "Pt");  add("pd",   "SQ",  2.05, "Pd")
add("ru",      "MH",  2.10, "Ru");  add("ir",   "MH",  2.10, "Ir")
add("au",      "LIN", 2.30, "Au")
# ---- toxicology ----
add("cd",      "TZ",  2.30, "Cd");  add("hg",   "LIN", 2.35, "Hg")
# ---- alkali / alkaline-earth ----
add("na",      "MH",  2.40, "Na");  add("k",    "MH",  2.70, "K")
add("li",      "MH",  2.10, "Li");  add("al",   "MH",  2.00, "Al")
add("sr",      "MH",  2.50, "Sr");  add("ba",   "MH",  2.70, "Ba")
# ---- early TMs ----
add("v",       "MH",  2.00, "V");   add("cr",   "MH",  2.00, "Cr")
add("ti",      "MH",  2.05, "Ti");  add("sc",   "MH",  2.10, "Sc")
add("y",       "MH",  2.25, "Y");   add("zr",   "MH",  2.15, "Zr")
add("nb",      "MH",  2.10, "Nb");  add("hf",   "MH",  2.15, "Hf")
add("ta",      "MH",  2.10, "Ta");  add("w",    "MH",  2.05, "W")
add("mo",      "MH",  2.05, "Mo")
# ---- 2nd/3rd-row TMs ----
add("rh",      "MH",  2.10, "Rh");  add("ag",   "LIN", 2.30, "Ag")
add("tc",      "MH",  2.10, "Tc");  add("re",   "MH",  2.10, "Re")
add("os",      "MH",  2.10, "Os")
# ---- post-transition / metalloids ----
add("ga",      "MH",  2.05, "Ga");  add("in",   "MH",  2.25, "In")
add("sn",      "MH",  2.30, "Sn");  add("sb",   "MH",  2.30, "Sb")
add("bi",      "MH",  2.40, "Bi");  add("tl",   "MH",  2.40, "Tl")
add("pb",      "MH",  2.45, "Pb")
# ---- s-block extras ----
add("rb",      "MH",  2.70, "Rb");  add("cs",   "MH",  2.80, "Cs")
add("ra",      "MH",  2.70, "Ra")
# ---- lanthanides ----
for tok, d in [("la",2.50),("ce",2.45),("pr",2.40),("nd",2.40),("pm",2.40),
               ("sm",2.35),("eu",2.40),("gd",2.35),("tb",2.35),("dy",2.35),
               ("ho",2.30),("er",2.30),("tm",2.30),("yb",2.30),("lu",2.25)]:
    add(tok, "MH", d, tok.capitalize())
# ---- metalloids with pharmacol. interest ----
add("se",      "MH",  2.30, "Se");  add("as",   "MH",  2.30, "As")
add("ge",      "MH",  2.10, "Ge");  add("b",    "MH",  2.05, "B")
add("be",      "MH",  2.05, "Be");  add("te",   "MH",  2.30, "Te")
add("po",      "MH",  2.40, "Po");  add("at",   "MH",  2.30, "At")
# ---- actinides ----
add("ac",      "MH",  2.50, "Ac");  add("th",   "MH",  2.45, "Th")
add("pa",      "MH",  2.40, "Pa");  add("u",    "MH",  2.40, "U")
add("np",      "MH",  2.35, "Np");  add("pu",   "MH",  2.40, "Pu")
add("am",      "MH",  2.40, "Am");  add("cm",   "MH",  2.35, "Cm")
add("bk",      "MH",  2.35, "Bk");  add("cf",   "MH",  2.30, "Cf")
add("es",      "MH",  2.30, "Es");  add("fm",   "MH",  2.30, "Fm")
# ---- oxidation-state variants ----
add("fe2",     "MH",  2.10, "Fe");  add("fe3",   "MH",  2.00, "Fe")
add("cu1",     "LIN", 2.15, "Cu");  add("cu2",   "SQ",  2.00, "Cu")
add("cu2_jt",  "JT",  2.00, "Cu")
add("mn2",     "MH",  2.20, "Mn");  add("mn3",   "MH",  2.00, "Mn")
add("mn3_jt",  "JT",  2.10, "Mn")
add("co2",     "MH",  2.05, "Co");  add("co3",   "MH",  2.00, "Co")
add("ni2",     "SQ",  1.95, "Ni");  add("ni3",   "SQ",  1.95, "Ni")
add("as3",     "MH",  2.30, "As");  add("as5",   "MH",  2.30, "As")
add("sb3",     "MH",  2.30, "Sb");  add("sb5",   "MH",  2.30, "Sb")
add("v4",      "MH",  2.00, "V");   add("v5",    "MH",  2.00, "V")
add("mo4",     "MH",  2.05, "Mo");  add("mo6",   "MH",  2.00, "Mo")
add("uo2",     "MH",  2.44, "U")
# ---- lanthanide chelate modes ----
add("gd_dtpa", "MH",  2.40, "Gd");  add("gd_dota","MH", 2.40, "Gd")
# ---- aqua modes ----
add("mg_aq",   "MH",  2.05, "Mg");  add("ca_aq",  "MH", 2.35, "Ca")
add("fe3_aq",  "MH",  2.00, "Fe");  add("mn2_aq", "MH", 2.20, "Mn")
add("co2_aq",  "MH",  2.05, "Co")

MH_DIRS = [(1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0),(0,0,-1)]
TZ_DIRS = [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
SQ_DIRS = [(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)]
LIN_DIRS= [(1,0,0),(-1,0,0)]

def run(cmd, timeout=180):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"

def first_model(path):
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

def score_only(rec, lig, mode):
    rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", mode,
                      "--receptor", rec, "--ligand", lig, "--score_only",
                      "--center_x", "0", "--center_y", "0", "--center_z", "0",
                      "--size_x", "16", "--size_y", "16", "--size_z", "16",
                      "--verbosity", "0", "--cpu", "4"])
    if rc != 0: return None
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[0-9.]+)", so)
    return float(m.group(1)) if m else None

def main():
    modes = sorted(M.keys())
    results = []
    for mode in modes:
        cls, d0, tok = M[mode]
        if cls == "TZ":   full = TZ_DIRS
        elif cls == "SQ": full = SQ_DIRS
        elif cls == "LIN":full = LIN_DIRS
        elif cls == "JT": full = SQ_DIRS      # 4 equatorial; axial weak
        else:             full = MH_DIRS

        enter = norm(full[-1])
        rec_donors = [tuple(c*d0 for c in norm(v)) for v in full[:-1]]
        d = os.path.join(SYN, mode); os.makedirs(d, exist_ok=True)
        rec = os.path.join(d, "rec.pdbqt")
        lig_ideal = os.path.join(d, "lig_ideal.pdbqt")
        lig_disp  = os.path.join(d, "lig_disp.pdbqt")
        write_receptor(rec, tok, rec_donors)
        write_ligand(lig_ideal, enter, d0)
        write_ligand(lig_disp, enter, d0 + 4.0)

        e_ideal = score_only(rec, lig_ideal, mode)
        e_far   = score_only(rec, lig_disp, mode)

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

        vout = os.path.join(d, "vina_out.pdbqt")
        rc2, so2, se2 = run([VINA, "--receptor", rec, "--ligand", lig_disp,
                             "--center_x", "0", "--center_y", "0", "--center_z", "0",
                             "--size_x", "12", "--size_y", "12", "--size_z", "12",
                             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                             "--out", vout])
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
        print(f"{mode:8s} {tok:3s} {cls:3s} d0={d0:.2f} | E_ideal={e_ideal if e_ideal is not None else 99:6.1f} "
              f"E_far={e_far if e_far is not None else 99:6.1f} | LKina d={d_lig if d_lig else 99:4.2f} "
              f"{'OK' if rc1==0 else 'FAIL'} | Vina d={vina_d if vina_d else 99:4.2f} "
              f"{'OK' if rc2==0 else 'FAIL'}")
        sys.stdout.flush()

    with open(os.path.join(BENCH, "metal_coverage_results_all.json"), "w") as f:
        json.dump(results, f, indent=2)

    n = len(results)
    n_ok   = sum(1 for r in results if r["lkina_rc"] == 0)
    n_vok  = sum(1 for r in results if r["vina_rc"] == 0)
    n_well = sum(1 for r in results if r["e_ideal"] is not None and r["e_far"] is not None
                 and r["e_ideal"] < r["e_far"] + 0.3)
    n_coord = sum(1 for r in results if r["coord"])
    print(f"\n=== SUMMARY: {n} modes | LKina-dock-OK {n_ok}/{n} | Vina-OK {n_vok}/{n} | "
          f"metal-well {n_well}/{n} | METAL_COORD reported {n_coord}/{n} ===")

if __name__ == "__main__":
    main()
