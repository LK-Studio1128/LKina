#!/usr/bin/env python3
"""Prepare hotspot-case systems (KRAS G12C 6OIM, CA II 3HS4) for LKina testing.

Outputs per system:
  <ID>_lig.pdb      co-crystallized ligand (heavy atoms, cleaned)
  <ID>_lig.pdbqt    ligand PDBQT (obabel, rigid redock style)
  <ID>_rec.pdb      receptor protein + catalytic metal (no waters/other HET)
  <ID>_rec.pdbqt    receptor PDBQT (obabel -xr, metal kept)
  <ID>_meta.json    box center/span from ligand, metal coords, ligand name
"""
import json
import subprocess
from pathlib import Path

BASE = Path(__file__).parent

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
         "HSD HSE HSP MSE".split())
METAL_ELEMS = {"ZN", "FE", "CU", "MG", "MN", "CA", "NI", "CO", "PT", "CD"}


def parse_pdb(path):
    atoms = []
    for ln in open(path, errors="ignore"):
        if ln.startswith(("ATOM", "HETATM")):
            atoms.append({
                "record": ln[:6].strip(),
                "name": ln[12:16].strip(),
                "resn": ln[17:20].strip(),
                "chain": ln[21].strip(),
                "resi": ln[22:26].strip(),
                "x": float(ln[30:38]), "y": float(ln[38:46]), "z": float(ln[46:54]),
                "elem": ln[76:78].strip(),
                "line": ln.rstrip(),
            })
    return atoms


def dist(a, b):
    return ((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2) ** 0.5


def prep(pdb_id, lig_resn, metal_resn):
    atoms = parse_pdb(BASE / f"{pdb_id}_raw.pdb")

    # ligand atoms: restrict to first occurrence (chain,resi) of lig_resn
    lig_all = [a for a in atoms if a["resn"] == lig_resn]
    keys = []
    for a in lig_all:
        k = (a["chain"], a["resi"])
        if k not in keys:
            keys.append(k)
    lig_key = keys[0]
    lig = [a for a in lig_all if (a["chain"], a["resi"]) == lig_key]
    print(f"  ligand copies found: {keys}; using {lig_key}")
    # metals
    metals = [a for a in atoms if a["resn"].strip() == metal_resn.strip()
              or (a["record"] == "HETATM" and a["elem"] in METAL_ELEMS
                  and a["resn"].strip() == metal_resn.strip())]
    if not metals:
        metals = [a for a in atoms if a["record"] == "HETATM"
                  and a["resn"].strip() == metal_resn.strip()]

    # protein: standard residues only, drop waters
    prot = [a for a in atoms if a["resn"] in AA]

    # write ligand pdb (heavy atoms)
    lp = BASE / f"{pdb_id}_lig.pdb"
    with open(lp, "w") as f:
        serial = 1
        for a in lig:
            if a["elem"] == "H":
                continue
            ln = a["line"]
            ln = f"{ln[:6]:6s}{serial:5d}{ln[11:]}" if len(ln) > 11 else ln
            f.write(ln + "\n")
            serial += 1
        f.write("END\n")

    # write receptor pdb: protein + metal
    rp = BASE / f"{pdb_id}_rec.pdb"
    with open(rp, "w") as f:
        serial = 1
        for a in prot:
            ln = a["line"]
            ln = f"{ln[:6]:6s}{serial:5d}{ln[11:]}"
            f.write(ln + "\n")
            serial += 1
        for a in metals:
            ln = a["line"]
            # force ZN resname for metal so obabel keeps element
            ln = f"HETATM{serial:5d} {a['name']:.>3s} {a['resn']:.>3s} {a['chain']}{a['resi']:>4s}" \
                 f"    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{1.00:6.2f}{20.00:6.2f}" \
                 f"          {a['elem']:>2s}  "
            f.write(ln + "\n")
            serial += 1
        f.write("END\n")

    # box from ligand center + span
    xs = [a["x"] for a in lig]; ys = [a["y"] for a in lig]; zs = [a["z"] for a in lig]
    center = [sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)]
    span = [max(xs)-min(xs), max(ys)-min(ys), max(zs)-min(zs)]

    meta = {
        "pdb": pdb_id, "lig": lig_resn, "metal": metal_resn,
        "center": center, "span": span,
        "metal_coords": [[m["x"], m["y"], m["z"]] for m in metals],
        "n_lig_atoms": len(lig), "n_prot_atoms": len(prot),
    }
    json.dump(meta, open(BASE / f"{pdb_id}_meta.json", "w"), indent=1)

    # convert: ligand pdbqt (rigid), receptor pdbqt (-xr keeps as receptor)
    subprocess.run(f'obabel "{lp}" -O "{BASE / f"{pdb_id}_lig.pdbqt"}" 2>/dev/null', shell=True)
    subprocess.run(f'obabel "{rp}" -O "{BASE / f"{pdb_id}_rec.pdbqt"}" -xr 2>/dev/null', shell=True)

    for check in ("lig.pdbqt", "rec.pdbqt"):
        p = BASE / f"{pdb_id}_{check}"
        n = sum(1 for _ in open(p)) if p.exists() else 0
        print(f"  {p.name}: {n} lines")
    print(f"  metal at {[(round(x,2),round(y,2),round(z,2)) for x,y,z in meta['metal_coords']]}")


if __name__ == "__main__":
    print("== 6OIM KRAS G12C (lig MOV, metal MG structural) ==")
    prep("6OIM", "MOV", "MG")
    print("== 3HS4 CA II (lig AZM, metal ZN catalytic) ==")
    prep("3HS4", "AZM", "ZN")
