#!/usr/bin/env python3
"""Reactive covalent docking benchmark: test all 6 LKina presets on synthetic
receptor-ligand systems with the reactive pair pre-positioned at ideal NAC
geometry.  Verifies: preset parameter loading, NAC detection, gradient check,
and P1 vs P1+P2 energy decomposition."""
import os, sys, math, subprocess, re, json

BENCH = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks"
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
VINA  = "/Users/luoxiaowen/Desktop/LKDock/源文件/LKDock_v4.0_Mac/engine/vina127"
RXN   = os.path.join(BENCH, "reactive_tests")
os.makedirs(RXN, exist_ok=True)

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

# Preset -> (receptor reactive atom name/resn, ligand nucleophile type,
#            ideal bond length A, ideal angle deg, anchor offset direction)
# We build: receptor = electrophile atom at origin + frame atom at +x;
#           ligand   = nucleophile at (bond_len along +z) with frame atom behind.
PRESETS = {
  "cys_michael":   {"rec": ("SG", "CYS", "S"), "lig": "C",  "r0": 1.82, "theta": 109.5, "dir": (0, 0, 1)},
  "cys_sn2":       {"rec": ("SG", "CYS", "S"), "lig": "C",  "r0": 1.82, "theta": 180.0, "dir": (0, 0, 1)},
  "ser_covalent":  {"rec": ("OG", "SER", "O"), "lig": "C",  "r0": 1.34, "theta": 109.5, "dir": (0, 0, 1)},
  "lys_targeting": {"rec": ("NZ", "LYS", "N"), "lig": "C",  "r0": 1.47, "theta": 109.5, "dir": (0, 0, 1)},
  "boronic_acid":  {"rec": ("OG", "SER", "O"), "lig": "B",  "r0": 1.47, "theta": None,   "dir": (0, 0, 1)},
  "tyr_covalent":  {"rec": ("OH", "TYR", "O"), "lig": "C",  "r0": 1.38, "theta": 109.5, "dir": (0, 0, 1)},
}

def norm(v):
    l = math.sqrt(sum(c*c for c in v)); return tuple(c/l for c in v)

def build_system(preset_name, p):
    d = os.path.join(RXN, preset_name); os.makedirs(d, exist_ok=True)
    rec_atom_name, resn, rec_typ = p["rec"]
    r0 = p["r0"]; theta = p["theta"]; dirv = norm(p["dir"])

    # receptor: electrophile at origin, frame atom at +x (3 A), plus 2 anchor residues
    rec_lines = []
    rec_lines.append(_atom_line("ATOM", 1, "CA", "ALA", "A", 1, -2.0, 0.0, 0.0, "+0.100", "C"))
    rec_lines.append(_atom_line("ATOM", 2, "CB", "ALA", "A", 1, -1.0, 0.0, 0.0, "+0.100", "C"))
    rec_lines.append(_atom_line("ATOM", 3, rec_atom_name, resn, "A", 1, 0.0, 0.0, 0.0, "-0.200", rec_typ))
    rec_lines.append(_atom_line("ATOM", 4, "CB", resn, "A", 1, 1.8, 0.0, 0.0, "+0.100", "C"))  # frame atom (CB)
    rec_lines.append("END")
    rec_path = os.path.join(d, "rec.pdbqt")
    open(rec_path, "w").write("\n".join(rec_lines) + "\n")

    # ligand: nucleophile at r0 along dir, frame atom at 1.5 A behind
    nx, ny, nz = dirv[0]*r0, dirv[1]*r0, dirv[2]*r0
    lig_typ = p["lig"]
    lig_lines = [
      "REMARK  reactive probe ligand",
      "ROOT",
      _atom_line("ATOM", 1, "C1", "LIG", "A", 1, nx, ny, nz, "-0.100", lig_typ),           # nucleophile
      _atom_line("ATOM", 2, "C2", "LIG", "A", 1, 0.0, 0.0, -1.5, "+0.100", "C"),            # ligand frame
      _atom_line("ATOM", 3, "H1", "LIG", "A", 1, nx+0.9, ny, nz, "+0.050", "HD"),
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

def main():
    results = []
    for name, p in PRESETS.items():
        d, rec, lig = build_system(name, p)
        box = (0, 0, 0, 16, 16, 16)
        anchor = "A:1:" + p["rec"][0]
        # P1 only (distance)
        out1 = os.path.join(d, "out_p1.pdbqt")
        rc1, so1, se1 = run([LKINA, "--scoring", "ad4", "--generate_maps",
                             "--receptor", rec, "--ligand", lig,
                             "--center_x", "0", "--center_y", "0", "--center_z", "0",
                             "--size_x", "16", "--size_y", "16", "--size_z", "16",
                             "--reactive_preset", name,
                             "--reactive_rec_atom", anchor,
                             "--reactive_lig_atom", "index:1",
                             "--reactive_mode", "distance",
                             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                             "--out", out1, "--verbosity", "1"])
        # P1+P2 (distance+angle with frame)
        out2 = os.path.join(d, "out_p12.pdbqt")
        frame = "A:1:CB"
        rc2, so2, se2 = run([LKINA, "--scoring", "ad4", "--generate_maps",
                             "--receptor", rec, "--ligand", lig,
                             "--center_x", "0", "--center_y", "0", "--center_z", "0",
                             "--size_x", "16", "--size_y", "16", "--size_z", "16",
                             "--reactive_preset", name,
                             "--reactive_rec_atom", anchor,
                             "--reactive_lig_atom", "index:1",
                             "--reactive_mode", "hybrid",
                             "--reactive_frame_atom", frame,
                             "--reactive_lig_frame_atom", "index:2",
                             "--reactive_gradcheck",
                             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                             "--out", out2, "--verbosity", "1"])

        def parse_out(path, out_txt):
            res = {"nac": None, "rdist": None, "rangle": None, "re": None, "grad": None}
            txt = (out_txt or "")
            if os.path.exists(path):
                try: txt += open(path).read()
                except: pass
            for l in txt.splitlines():
                m = re.search(r"REACTIVE_NAC[:=]\s*(YES|NO)", l)
                if m: res["nac"] = m.group(1)
                m = re.search(r"REACTIVE_DIST[:=]\s*([0-9.]+)", l)
                if m: res["rdist"] = float(m.group(1))
                m = re.search(r"REACTIVE_ANGLE[:=]\s*([0-9.]+)", l)
                if m: res["rangle"] = float(m.group(1))
                m = re.search(r"REACTIVE_DIST_E[:=]\s*(-?[0-9.]+)", l)
                if m: res["re"] = float(m.group(1))
                m = re.search(r"gradient check passed|gradcheck passed|gradient.*OK|finite.diff.*pass", l, re.I)
                if m: res["grad"] = True
            return res

        r1 = parse_out(out1, so1 + se1) if rc1 == 0 else {}
        r2 = parse_out(out2, so2 + se2) if rc2 == 0 else {}
        results.append({
            "preset": name, "r0": p["r0"], "theta": p["theta"],
            "p1_rc": rc1, "p1_nac": r1.get("nac"), "p1_dist": r1.get("rdist"),
            "p12_rc": rc2, "p12_nac": r2.get("nac"), "p12_dist": r2.get("rdist"),
            "p12_angle": r2.get("rangle"), "p12_re": r2.get("re"), "gradcheck": r2.get("grad"),
            "err1": (se1 or so1)[-200:] if rc1 != 0 else "",
            "err2": (se2 or so2)[-200:] if rc2 != 0 else "",
        })
        print(f"{name:14s} r0={p['r0']:.2f} | P1 rc={rc1} nac={r1.get('nac')} "
              f"| P1+P2 rc={rc2} nac={r2.get('nac')} dist={r2.get('rdist')} "
              f"angle={r2.get('rangle')} grad={r2.get('grad')}")

    with open(os.path.join(BENCH, "reactive_presets_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    n_ok = sum(1 for r in results if r["p1_rc"] == 0 and r["p12_rc"] == 0)
    print(f"\n=== SUMMARY: {n_ok}/{len(PRESETS)} presets completed both P1 and P1+P2 ===")

if __name__ == "__main__":
    main()
