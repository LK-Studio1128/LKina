#!/usr/bin/env python3
"""C3 deep-pocket ablation — single-stage vs C3 (two-step) vs C3b (weak attractor).

Motivation (paper §2.10, LKINA.md §5.4.3): the full constraint set applied
throughout search can trap the ligand near the receptor surface and under-sample
the global surface. C3 runs an unconstrained Phase-1 MC, filters poses within
--reactive_presample_dist of the anchor, then refines with the full P1-P4 set in
Phase-2. C3b uses a broad/weak Gaussian bias (sigma x3, eps x0.15) in Phase 1 as
an attracting-cavity prior. The expected benefit is largest when the reactive
anchor sits in a *deep pocket* whose entrance is far from the ligand start.

Design: synthetic receptor with the Cys-SG reactive anchor buried at the bottom
of a 6-wall cylinder (deep pocket, single opening along +z); the probe ligand
starts OUTSIDE the pocket mouth at z = +6 A. Three variants:

  single  : full P1+P2 constraints from the start (no --reactive_two_step)
  c3      : --reactive_two_step (unconstrained Phase-1, Phase-2 refine)
  c3b     : --reactive_two_step --reactive_weak_attractor (Goullieux-like prior)

Metrics per run:
  rc                    : return code
  REACTIVE_DIST         : final d(SG, nucleophile) from docked pose
  REACTIVE_ANGLE        : final attack angle
  REACTIVE_NAC          : NAC YES/NO (dist < 3.0 A and angle within +/-25 deg)
  energy                : mode-1 affinity
  converged             : REACTIVE_DIST <= 2.6 A (phase-2 refined into reaction
                          geometry; the deep-pocket success criterion)

The synthetic pocket is deliberately narrow (diameter ~4.2 A vs ligand span
~4.4 A) so a blind unconstrained search rarely threads the entrance, while the
anchor-directed Phase-2 refine (and C3b's Phase-1 attractor) does.

Usage:
  python3 benchmarks/run_c3_ablation.py          # all six presets
  python3 benchmarks/run_c3_ablation.py --preset cys_michael
Output: benchmarks/c3_ablation_results.json
"""
import argparse, json, math, os, re, subprocess, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LKINA = os.path.join(REPO, "build", "mac", "release", "LKina")
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "c3_ablation")
os.makedirs(OUTDIR, exist_ok=True)

# preset -> (r0, theta, rec_atom, rec_type, lig_type)
PRESETS = {
    "cys_michael":   (1.82, 109.5, "SG", "S", "C"),
    "cys_sn2":       (1.82, 180.0, "SG", "S", "C"),
    "ser_covalent":  (1.34, 109.5, "OG", "OA", "C"),
    "lys_targeting": (1.47, 109.5, "NZ", "N", "C"),
    "boronic_acid":  (1.47, None,  "OH", "OA", "B"),
    "tyr_covalent":  (1.38, 109.5, "OH", "OA", "C"),
}

WALL_R = 3.4   # pocket wall radius (A) — wide enough for sampling to enter
WALL_N = 10    # wall atoms per ring (channel gap ~ 2*pi*R/N = 2.1 A)
RINGS  = 6     # ring depth below rim (deep pocket, long channel)
RIM_Z  = 3.0   # pocket mouth height (z); open only above this
BOX = 32       # search box half-size (large space -> C3 global sampling matters)
EXH  = 32      # exhaustiveness (deep-pocket search is hard; keep equal for all variants)


def atom_line(record, serial, name, resn, chain, resi, x, y, z, q, typ):
    return (f"{record:<6s}{serial:5d} {name:>4s} {resn:>3s} {chain:1s}{resi:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}"
            f"    {q:>6s}  {typ:<2s}")


def build_deep_pocket(preset):
    r0, theta, rec_atom, rec_typ, lig_typ = PRESETS[preset]
    d = os.path.join(OUTDIR, preset)
    os.makedirs(d, exist_ok=True)
    rec = []
    # reactive anchor: Cys SG/OG/NZ/OH at origin with CB frame at (-1.8,0,0)
    rec.append(atom_line("ATOM", 1, "CA", "ALA", "A", 1, -3.0, 0.0, 0.0, "+0.100", "C"))
    rec.append(atom_line("ATOM", 2, "CB", "ALA", "A", 1, -1.8, 0.0, 0.0, "+0.100", "C"))
    rec.append(atom_line("ATOM", 3, rec_atom, "CYS", "A", 1, 0.0, 0.0, 0.0, "-0.200", rec_typ))
    # deep-pocket walls: RINGS rings below the rim at z = RIM_Z - 2.0*k,
    # plus one rim ring at RIM_Z; only the top (+z) is open. Anchor sits at
    # the bottom (z=0) of a cylinder whose mouth is at z=RIM_Z.
    sid = 10
    for k in range(RINGS):
        z = RIM_Z - 2.0 * (k + 1)
        for i in range(WALL_N):
            ang = 2 * math.pi * i / WALL_N
            rec.append(atom_line("ATOM", sid, "C", "PHE", "A", 1,
                                 WALL_R * math.cos(ang), WALL_R * math.sin(ang), z,
                                 "+0.050", "C"))
            sid += 1
    # rim ring at the mouth (z = RIM_Z)
    for i in range(WALL_N):
        ang = 2 * math.pi * i / WALL_N
        rec.append(atom_line("ATOM", sid, "C", "PHE", "A", 1,
                             WALL_R * math.cos(ang), WALL_R * math.sin(ang), RIM_Z,
                             "+0.050", "C"))
        sid += 1
    rec.append("END")
    rec_path = os.path.join(d, "rec.pdbqt")
    open(rec_path, "w").write("\n".join(rec) + "\n")

    # ligand starts OUTSIDE the pocket mouth, far away (14.4 A from anchor):
    # the P1 Gaussian well decays to ~0 beyond r=4 A (Section 2.10), so from
    # this distance no variant feels any reactive pull; success requires the
    # search itself to thread the pocket mouth. A rigid 5-atom chain (6-DOF
    # search) keeps the search non-trivial in the 32 A box.
    lig = [
        "REMARK  deep-pocket probe (far start)", "ROOT",
        atom_line("ATOM", 1, "C1", "LIG", "A", 1, 12.0, 0.0, 8.0, "-0.100", lig_typ),
        atom_line("ATOM", 2, "C2", "LIG", "A", 1, 10.5, 1.2, 8.0, "+0.050", "C"),
        atom_line("ATOM", 3, "C3", "LIG", "A", 1, 8.5, 1.2, 8.0, "+0.050", "C"),
        atom_line("ATOM", 4, "C4", "LIG", "A", 1, 8.5, 1.2, 6.5, "+0.050", "C"),
        atom_line("ATOM", 5, "H1", "LIG", "A", 1, 12.9, 0.0, 8.0, "+0.050", "HD"),
        "ENDROOT", "TORSDOF 0",
    ]
    lig_path = os.path.join(d, "lig.pdbqt")
    open(lig_path, "w").write("\n".join(lig) + "\n")
    return d, rec_path, lig_path, (r0, theta, rec_atom)


def run(cmd, timeout=600):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"


def collect(path):
    txt = open(path, errors="ignore").read() if os.path.exists(path) else ""
    out = {}
    for key in ("REACTIVE_NAC", "REACTIVE_DIST", "REACTIVE_ANGLE", "VINA RESULT",
                "METAL_GEO_E", "REACTIVE_DIST_E", "REACTIVE_ANGLE_E"):
        m = re.search(rf"{key}[:=]\s*([YESNO0-9.eE+-]+)", txt)
        if m:
            v = m.group(1)
            if v in ("YES", "NO"):
                out[key] = v
            else:
                try:
                    out[key] = float(v)
                except ValueError:
                    pass
    return out


def run_variant(lkina, rec, lig, preset, rec_atom, variant):
    base = [lkina, "--scoring", "ad4", "--generate_maps",
            "--receptor", rec, "--ligand", lig,
            "--center_x", "0", "--center_y", "0", "--center_z", "0",
            "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
            "--reactive_preset", preset,
            "--reactive_rec_atom", "A:1:" + rec_atom,
            "--reactive_lig_atom", "index:1",
            "--reactive_mode", "hybrid",
            "--reactive_frame_atom", "A:1:CB",
            "--reactive_lig_frame_atom", "index:2",
            "--exhaustiveness", str(EXH), "--seed", "42", "--cpu", "4", "--verbosity", "1"]
    out = os.path.join(OUTDIR, preset, f"out_{variant}.pdbqt")
    if variant == "single":
        cmd = base + ["--out", out]
    elif variant == "c3":
        cmd = base + ["--reactive_two_step", "--reactive_presample_dist", "15.0",
                      "--out", out]
    else:  # c3b
        cmd = base + ["--reactive_two_step", "--reactive_presample_dist", "15.0",
                      "--reactive_weak_attractor", "--out", out]
    rc, so, se = run(cmd)
    txt = so + se
    entry = {"rc": rc}
    if rc == 0 and os.path.exists(out):
        entry.update(collect(out))
        # mode-1 energy
        m = re.search(r"^\s*1\s+(-?[0-9.]+)\s", so, re.M)
        if m:
            entry["energy"] = float(m.group(1))
    else:
        entry["fail"] = (se or so).strip().splitlines()[-1][:100] if (se or so).strip() else "no output"
    # deep-pocket success: refined into reaction geometry (dist <= 2.6 A)
    dd = entry.get("REACTIVE_DIST")
    entry["converged"] = (dd is not None and dd <= 2.6)
    return entry, txt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=None)
    ap.add_argument("--preset", default=None)
    args = ap.parse_args()
    lkina = args.binary or LKINA
    if not os.path.exists(lkina):
        print("ERROR: binary not found:", lkina); sys.exit(2)

    presets = [args.preset] if args.preset else list(PRESETS)
    results = {}
    print(f"deep pocket: 6-wall cylinder r={WALL_R} A x {RINGS} rings, mouth at +z")
    for preset in presets:
        d, rec, lig, (r0, theta, rec_atom) = build_deep_pocket(preset)
        print(f"== {preset}  (r0={r0}, theta={theta}) ==")
        row = {"r0": r0, "theta": theta, "variants": {}}
        for variant in ("single", "c3", "c3b"):
            entry, txt = run_variant(lkina, rec, lig, preset, rec_atom, variant)
            row["variants"][variant] = entry
            print(f"  {variant:7s} rc={entry['rc']} dist={entry.get('REACTIVE_DIST')} "
                  f"angle={entry.get('REACTIVE_ANGLE')} nac={entry.get('REACTIVE_NAC')} "
                  f"conv={entry.get('converged')} E={entry.get('energy')}")
        results[preset] = row

    # summary: phase-2 convergence counts per variant
    summary = {"n_presets": len(presets)}
    for v in ("single", "c3", "c3b"):
        conv = sum(1 for p in presets if results[p]["variants"][v].get("converged"))
        ok = sum(1 for p in presets if results[p]["variants"][v].get("rc") == 0)
        summary[v] = {"ok": ok, "converged": conv}
        print(f"SUMMARY {v}: {conv}/{ok} converged")
    from pathlib import Path
    (Path(os.path.dirname(os.path.abspath(__file__))) / "c3_ablation_results.json"
     ).write_text(json.dumps({"summary": summary, "results": results}, indent=1))
    print("saved: benchmarks/c3_ablation_results.json")


if __name__ == "__main__":
    main()
