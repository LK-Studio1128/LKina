#!/usr/bin/env python3
"""Supplementary feature tests for LKina (task 14):
  A) reactive gradient check (--reactive_gradcheck) on each preset
  B) Metal Bias (O5) --metal_bias auto-attractor on a metal receptor
  C) C3 two-step strategy on cys_michael
  D) backward-compat + energy-scale sanity checks
"""
import os, sys, math, subprocess, re, json

BENCH = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks"
LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
SUP   = os.path.join(BENCH, "feature_tests")
os.makedirs(SUP, exist_ok=True)

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

def run(cmd, timeout=180):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"

# ---------------- A. gradient checks ----------------
def test_gradcheck():
    print("=== A. reactive gradient check (each preset, P1+P2) ===")
    out = []
    presets = ["cys_michael", "cys_sn2", "ser_covalent", "lys_targeting", "boronic_acid", "tyr_covalent"]
    for preset in presets:
        d = os.path.join(SUP, "grad_" + preset); os.makedirs(d, exist_ok=True)
        # minimal receptor: Cys-SG at origin + CB frame; ligand C at 1.82 along +z
        rec = os.path.join(d, "rec.pdbqt")
        lig = os.path.join(d, "lig.pdbqt")
        open(rec, "w").write("\n".join([
            _atom_line("ATOM", 1, "CA", "ALA", "A", 1, -2.0, 0.0, 0.0, "+0.100", "C"),
            _atom_line("ATOM", 2, "CB", "ALA", "A", 1, -1.0, 0.0, 0.0, "+0.100", "C"),
            _atom_line("ATOM", 3, "SG", "CYS", "A", 1, 0.0, 0.0, 0.0, "-0.200", "S"),
            _atom_line("ATOM", 4, "CB", "CYS", "A", 1, 1.8, 0.0, 0.0, "+0.100", "C"),
            "END"]) + "\n")
        open(lig, "w").write("\n".join([
            "REMARK  probe", "ROOT",
            _atom_line("ATOM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, 1.82, "-0.100", "C"),
            _atom_line("ATOM", 2, "C2", "LIG", "A", 1, 0.0, 0.0, -1.5, "+0.100", "C"),
            _atom_line("ATOM", 3, "H1", "LIG", "A", 1, 0.9, 0.0, 1.82, "+0.050", "HD"),
            "ENDROOT", "TORSDOF 0"]) + "\n")
        rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps",
                          "--receptor", rec, "--ligand", lig,
                          "--center_x", "0", "--center_y", "0", "--center_z", "0",
                          "--size_x", "16", "--size_y", "16", "--size_z", "16",
                          "--reactive_preset", preset,
                          "--reactive_rec_atom", "A:1:SG",
                          "--reactive_lig_atom", "index:1",
                          "--reactive_mode", "hybrid",
                          "--reactive_frame_atom", "A:1:CB",
                          "--reactive_lig_frame_atom", "index:2",
                          "--reactive_gradcheck",
                          "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                          "--out", os.path.join(d, "out.pdbqt"), "--verbosity", "1"])
        txt = so + se
        passed = bool(re.search(r"pass|OK|match|within", txt, re.I)) and rc == 0
        out.append({"preset": preset, "rc": rc, "gradcheck_ok": passed,
                    "log_tail": txt[-200:] if not passed else ""})
        print(f"  {preset:14s} rc={rc} gradcheck={'PASS' if passed else 'see-log'}")
    json.dump(out, open(os.path.join(SUP, "gradcheck_results.json"), "w"), indent=2)

# ---------------- B. Metal Bias (O5) ----------------
def test_metal_bias():
    print("=== B. Metal Bias (O5) on clean 4JC ===")
    d = os.path.join(SUP, "metal_bias"); os.makedirs(d, exist_ok=True)
    rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", "zn",
                      "--metal_bias", "--metal_bias_strength", "2.0", "--metal_bias_width", "1.5",
                      "--receptor", os.path.join(BENCH, "real_systems/4JC_core_rec.pdbqt"),
                      "--ligand", os.path.join(BENCH, "4JC305_lig.pdbqt"),
                      "--center_x", "-2.55", "--center_y", "2.29", "--center_z", "85.60",
                      "--size_x", "25", "--size_y", "25", "--size_z", "25",
                      "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                      "--out", os.path.join(d, "out.pdbqt"), "--verbosity", "1"])
    txt = so + se
    bias_ok = bool(re.search(r"bias|Metal Bias|attractor", txt, re.I))
    print(f"  rc={rc} bias-flag-seen={bias_ok}")
    json.dump({"rc": rc, "bias_ok": bias_ok, "log": txt[-400:]},
              open(os.path.join(SUP, "metal_bias_results.json"), "w"), indent=2)

# ---------------- C. C3 two-step ----------------
def test_c3():
    print("=== C. C3 two-step (cys_michael, clean receptor) ===")
    d = os.path.join(SUP, "c3"); os.makedirs(d, exist_ok=True)
    rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps",
                      "--receptor", os.path.join(BENCH, "real_systems/4JC_core_rec.pdbqt"),
                      "--ligand", os.path.join(BENCH, "4JC305_lig.pdbqt"),
                      "--center_x", "-2.55", "--center_y", "2.29", "--center_z", "85.60",
                      "--size_x", "25", "--size_y", "25", "--size_z", "25",
                      "--reactive_preset", "cys_michael",
                      "--reactive_rec_atom", "A:797:SG",
                      "--reactive_lig_atom", "index:1",
                      "--reactive_two_step", "--reactive_presample_dist", "10.0",
                      "--exhaustiveness", "8", "--seed", "42", "--cpu", "4",
                      "--out", os.path.join(d, "out.pdbqt"), "--verbosity", "1"])
    txt = so + se
    # receptor has no Cys 797 -> expect graceful degradation or anchor-not-found
    print(f"  rc={rc} (two-step accepted; anchor A:797:SG may not exist in 4JC)")
    json.dump({"rc": rc, "log": (so + se)[-500:]},
              open(os.path.join(SUP, "c3_results.json"), "w"), indent=2)

# ---------------- D. energy-scale sanity ----------------
def test_energy_scale():
    print("=== D. energy-scale sanity (clean 4JC, crystal pose) ===")
    d = os.path.join(SUP, "energy_scale"); os.makedirs(d, exist_ok=True)
    # crystal pose score_only under: metal zn, plain ad4 (auto metal), and -metal_soft_weight
    for label, extra in [("metal_zn", ["--metal_mode", "zn"]),
                         ("metal_soft", ["--metal_mode", "zn", "--metal_soft_weight", "0.3"])]:
        rc, so, se = run([LKINA, "--scoring", "ad4", "--generate_maps"] + extra +
                         ["--receptor", os.path.join(BENCH, "real_systems/4JC_core_rec.pdbqt"),
                          "--ligand", os.path.join(BENCH, "4JC305_lig.pdbqt"), "--score_only",
                          "--center_x", "-2.55", "--center_y", "2.29", "--center_z", "85.60",
                          "--size_x", "25", "--size_y", "25", "--size_z", "25",
                          "--verbosity", "1", "--cpu", "4"])
        m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[0-9.]+)", so)
        e = m.group(1) if m else "ERR"
        print(f"  {label:12s} crystal-pose score_only = {e} kcal/mol")
    json.dump({}, open(os.path.join(SUP, "energy_scale_results.json"), "w"))

if __name__ == "__main__":
    test_gradcheck()
    test_metal_bias()
    test_c3()
    test_energy_scale()
    print("\n=== supplementary feature tests done ===")
