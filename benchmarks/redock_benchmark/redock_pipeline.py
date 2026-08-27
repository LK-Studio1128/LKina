#!/usr/bin/env python3
"""
LKina large-benchmark redocking pipeline.

Protocol (follows Santos-Martins et al. 2014 AutoDock4Zn redocking):
  1. Download PDB entries containing the target metal (ZN/FE/CU) with a bound ligand.
  2. Extract receptor (protein + metal) and co-crystallized ligand.
  3. Prepare PDBQT: receptor via obabel, ligand via obabel rigid re-dock.
  4. Redock ligand with:
       a) LKina AD4 + metal_mode (tz for Zn, mh for Fe3+, jt for Cu2+)
       b) LKina AD4 standard (no metal mode) -- proxy for "standard AD4"
       c) Vina 1.2.7 on identical LKina-generated AD4 maps (AD4 scoring)
     Note: true AutoDock4Zn parameters live inside LKina's zn mode.
  5. Score top pose vs crystal pose: top-1 RMSD <= 2.0 A = success.
"""
import json
import math
import os
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).parent
PDB_DIR = BASE / "pdb"
PREP_DIR = BASE / "prepared"
DOCK_DIR = BASE / "docking"
RES_DIR = BASE / "results"
for d in (PDB_DIR, PREP_DIR, DOCK_DIR, RES_DIR):
    d.mkdir(exist_ok=True)

LKEINA = "/tmp/lkina_release/LKina-macos-arm64"
VINA127 = "/Volumes/LK/软件资料/软件开发/分子对接全流程工具-PLIP版/vina"

# metal -> (metal_mode_token, ideal coordination element letter) per LKina CLI
METAL_MODES = {"ZN": "zn", "FE": "fe3", "CU": "cu2_jt"}


def sh(cmd, cwd=None, timeout=600):
    return subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True, timeout=timeout)


# ---------------------------------------------------------------- RCSB search

def rcsb_search(metal_code, max_results):
    """Entries with the given HET metal + at least one non-polymer organic ligand."""
    q = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                            "operator": "exact_match", "value": metal_code}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "less_or_equal", "value": 2.5}},
        ]},
        "request_options": {"paginate": {"start": 0, "rows": min(max_results, 500)},
                            "results_verbosity": "compact"},
        "return_type": "entry",
    }
    import urllib.request
    req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",
                                 data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    return d.get("identifier_set") or d.get("result_set") or [], int(d.get("total_count") or 0)


# ---------------------------------------------------------------- PDB parsing

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
         "HSD HSE HSP MSE SEP TPO PTR CSO CME OCS CAS MHO SEC PYL".split())


def parse_pdb(path):
    """Return atoms: list of dicts from ATOM/HETATM records."""
    atoms = []
    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                rec = line[:6].strip()
                name = line[12:16].strip()
                resn = line[17:20].strip()
                chain = line[21:22].strip()
                resi = line[22:26].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                elem = line[76:78].strip().upper() or "".join(c for c in name if c.isalpha())[:1]
                atoms.append({"rec": rec, "name": name, "resn": resn, "chain": chain,
                              "resi": resi, "xyz": (x, y, z), "elem": elem, "line": line})
    return atoms


def dist(a, b):
    return math.dist(a["xyz"], b["xyz"])


def select_ligand(atoms, metal_id):
    """Pick the largest HETATM residue within 6 A of any metal ion of matching elem."""
    metals = [a for a in atoms if a["resn"] == metal_id]
    if not metals:
        return None
    het_res = defaultdict(list)
    for a in atoms:
        if a["resn"] not in AA and a["resn"] != metal_id and len(a["name"]) > 0 and a["name"] != "O":
            # skip waters / ions
            if a["resn"] in ("HOH", "DOD", "WAT", "NA", "CL", "K", "MG", "CA", "SO4", "PO4", "GOL", "EDO", "ACT"):
                continue
            het_res[(a["chain"], a["resi"])].append(a)
    best, best_n = None, 0
    for key, group in het_res.items():
        elems = {a["elem"] for a in group}
        if not elems <= {"C", "N", "O", "S", "P", "F", "Cl", "BR", "I", "H"} | {"CL", "BR", "B", "SI", "SE"}:
            continue
        n_heavy = sum(1 for a in group if a["elem"] != "H")
        if n_heavy < 6:  # skip tiny fragments
            continue
        dmin = min(dist(l, m) for l in group for m in metals)
        if dmin <= 6.0 and n_heavy > best_n:
            best, best_n = group, n_heavy
    return best


def clean_pdb_write(lines, path, end="END"):
    body = "\n".join(ln.rstrip() for ln in lines if ln.strip())
    path.write_text(body + f"\n{end}\n")


def build_files(pdb_id, metal_code):
    """Write ligand pdb / receptor-with-metal pdb for one entry."""
    src = PDB_DIR / f"{pdb_id}.pdb"
    atoms = parse_pdb(src)
    lig = select_ligand(atoms, metal_code)
    if lig is None:
        return None
    metals = [a for a in atoms if a["resn"] == metal_code]
    lig_keys = {(a["rec"], a["chain"], a["resi"]) for a in lig}

    protein_metal = []
    seen = set()
    for a in atoms:
        key = (a["rec"], a["chain"], a["resi"])
        is_prot = a["resn"] in AA or a["rec"] == "ATOM"
        if is_prot or a["resn"] == metal_code or key in lig_keys:
            k2 = (a["serial"] if False else id(a))
            pass
    # simple: keep AA residues + chosen metals only for receptor; exclude ligand residue & waters
    lig_chain, lig_resi = lig[0]["chain"], lig[0]["resi"]
    prot_lines, metal_lines = [], []
    seen_res = set()
    for a in atoms:
        if a["resn"] in ("HOH", "DOD", "WAT"):
            continue
        if a["chain"] == lig_chain and a["resi"] == lig_resi and a["resn"] not in AA:
            continue
        if a["resn"] == metal_code:
            metal_lines.append(a["line"])
            continue
        if a["resn"] in AA:
            kk = (a["chain"], a["resi"])
            if kk in seen_res:
                continue
            seen_res.add(kk)
    # simpler approach: write full filtered atom lines
    lig_ids = {(a["chain"], a["resi"], a["resn"]) for a in lig}
    rec_lines = []
    seen_atoms = set()
    for a in atoms:
        lid = (a["chain"], a["resi"], a["resn"])
        if lid in lig_ids:
            continue
        if a["resn"] in ("HOH", "DOD", "WAT", "SO4", "PO4", "GOL", "EDO", "ACT", "NO3", "CL", "NA"):
            continue
        if a["resn"] in AA or a["rec"] == "ATOM":
            rec_lines.append(a["line"])
        elif a["resn"] == metal_code:
            rec_lines.append(a["line"])
    # deduplicate alternate conformations: keep first altloc per (chain,resi,name)
    out_recp, seen = [], set()
    for ln in rec_lines:
        altloc = ln[16]
        nm = (ln[21], ln[22:26], ln[12:16], ln[17:20])
        if altloc not in (" ", "A"):
            continue
        if nm in seen:
            continue
        seen.add(nm)
        out_recp.append(ln)
    out_lig, seenl = [], set()
    for ln in [a["line"] for a in lig]:
        altloc = ln[16]
        nm = (ln[12:16])
        if altloc not in (" ", "A"):
            continue
        if nm in seenl:
            continue
        seenl.add(nm)
        out_lig.append(ln)

    tag = f"{pdb_id}_{metal_code}"
    rp = PREP_DIR / f"{tag}_rec.pdb"
    lp = PREP_DIR / f"{tag}_lig.pdb"
    clean_pdb_write(out_lig, lp, end="TER\nEND")
    clean_pdb_write(out_recp, rp, end="END")
    # reference ligand coords (crystal pose)
    ref = RES_DIR / f"{tag}_ref.xyz"
    ref.write_text(str(len(out_lig)) + "\n\n")
    # parse to xyz list
    xyz = []
    for ln in out_lig:
        xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    tag_json = RES_DIR / f"{tag}_meta.json"
    met_coord = [(round(m["xyz"][0], 2), round(m["xyz"][1], 2), round(m["xyz"][2], 2)) for m in metals]
    json.dump({"pdb": pdb_id, "metal": metal_code, "metals": met_coord,
               "lig_center": [round(sum(p[i] for p in xyz) / len(xyz), 3) for i in range(3)],
               "span": [round(max(p[i] for p in xyz) - min(p[i] for p in xyz), 1) for i in range(3)],
               "lig_name": lig[0]["resn"]}, open(tag_json, "w"))
    return lp, rp, xyz


# ---------------------------------------------------------------- preparation

def prep_pdbqt(tag):
    lp = PREP_DIR / f"{tag}_lig.pdb"
    rp = PREP_DIR / f"{tag}_rec.pdb"
    lig_qt = PREP_DIR / f"{tag}_lig.pdbqt"
    rec_qt = PREP_DIR / f"{tag}_rec.pdbqt"
    if not lig_qt.exists() or lig_qt.stat().st_size == 0:
        r = sh(f'obabel "{lp}" -O "{lig_qt}" 2>/dev/null')
        if not lig_qt.exists() or lig_qt.stat().st_size == 0:
            return None
        txt = lig_qt.read_text()
        if "ROOT" not in txt:
            return None
        if txt.count("ROOT") > 1:  # disconnected multi-fragment ligand (e.g. pyrophosphate)
            return None
    if not rec_qt.exists() or rec_qt.stat().st_size == 0:
        r = sh(f'obabel "{rp}" -O "{rec_qt}" -xr 2>/dev/null')
        if not rec_qt.exists() or rec_qt.stat().st_size == 0:
            return None
        # obabel writes receptor as many ROOT/ENDROOT blocks => that's fine for rigid receptor in vina?
        # Vina requires no ROOT/TORSDOF for receptor. Clean them:
        lines = []
        for ln in rec_qt.read_text().splitlines():
            s = ln.strip()
            if s in ("ROOT", "ENDROOT", "TORSDOF 0", "BRANCH", "ENDBRANCH"):
                continue
            lines.append(ln)
        rec_qt.write_text("\n".join(lines) + "\n")
    return rec_qt, lig_qt


def rmsd_docked(docked_pdbqt, ref_xyz):
    """RMSD between first MODEL heavy atoms and reference order-matched by element sequence."""
    models = open(docked_pdbqt, errors="ignore").read().split("MODEL")
    first = None
    if len(models) > 1:
        first = models[1]
    else:
        first = open(docked_pdbqt, errors="ignore").read()
    dock = []
    for ln in first.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            name = ln[12:16].strip()
            elem = ln[77:79].strip().upper() or "".join(c for c in name if c.isalpha())[:1]
            if elem in ("H", "HD") or elem.startswith("H"):
                continue
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            adtype = ln[77:79].strip()
            mapelem = {"A": "C", "OA": "O", "NA": "N", "N": "N", "SA": "S", "S": "S",
                       "F": "F", "Cl": "Cl", "CL": "Cl", "Br": "Br", "BR": "Br", "I": "I",
                       "OS": "O", "NS": "N", "P": "P"}.get(adtype, elem)
            dock.append((mapelem, xyz))
    return dock


_AD_ELEM = {"A": "C", "OA": "O", "NA": "N", "N": "N", "NS": "N", "SA": "S", "S": "S",
            "OS": "O", "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P",
            "SI": "Si", "B": "B", "SE": "Se"}


def _heavy_atoms_pdbqt(path):
    """Ordered heavy-atom [(elem,(x,y,z))] from a ligand PDBQT (input or docked)."""
    out = []
    for ln in open(path, errors="ignore"):
        if ln.startswith(("ATOM", "HETATM")):
            adt = ln[77:79].strip().upper()
            name = ln[12:16].strip()
            elem = _AD_ELEM.get(adt)
            if elem is None:
                elem = "".join(c for c in name if c.isalpha())[:1].upper()
            if elem.startswith("H"):
                continue
            out.append((elem, (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))))
        elif ln.strip() in ("ENDROOT", "TORSDOF 0") :
            pass
    return out


def rmsd_best_effort(tag, docked_pdbqt):
    """Direct (no-superposition) RMSD, pairing heavy atoms in file order.

    Valid because the docked pose preserves the input ligand PDBQT atom order,
    so the k-th heavy atom of the pose corresponds to the k-th heavy atom of
    the co-crystallized ligand converted by the same obabel step."""
    try:
        ref = _heavy_atoms_pdbqt(PREP_DIR / f"{tag}_lig.pdbqt")
    except FileNotFoundError:
        return None
    txt = open(docked_pdbqt, errors="ignore").read()
    models = txt.split("MODEL")
    seg = models[1] if len(models) > 1 else txt
    tmp = Path("/tmp") / f"_pose_{tag}.pdbqt"
    tmp.write_text(seg if len(models) > 1 else seg.replace("ROOT", "", 1))
    dock = _heavy_atoms_pdbqt(tmp)
    tmp.unlink(missing_ok=True)
    if not ref or not dock or len(ref) != len(dock):
        return None
    if {e for e, _ in ref} != {e for e, _ in dock}:
        return None
    s = sum(math.dist(r, d) ** 2 for (_, r), (_, d) in zip(ref, dock)) / len(ref)
    return round(math.sqrt(s), 2)


# ---------------------------------------------------------------- docking

BOX_MARGIN = 6.0


def run_one(tag):
    prepped = prep_pdbqt(tag)
    if prepped is None:
        return {"id": tag, "status": "prep_fail"}
    rec_qt, lig_qt = prepped
    meta = json.load(open(RES_DIR / f"{tag}_meta.json"))
    cx, cy, cz = meta["lig_center"]
    sx, sy, sz = [max(s + BOX_MARGIN * 2, 10) for s in meta["span"]]
    metal = meta["metal"]
    mode = METAL_MODES[metal]
    res = {"id": tag, "status": "ok", "lig": meta["lig_name"]}

    box = ["--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", f"{sx:.0f}", "--size_y", f"{sy:.0f}", "--size_z", f"{sz:.0f}"]

    def lkina_cmd(extra, out):
        return ([LKEINA, "--scoring", "ad4", "--generate_maps"] + extra +
                ["--receptor", str(rec_qt), "--ligand", str(lig_qt)] + box +
                ["--out", out, "--seed", "42", "--exhaustiveness", "8",
                 "--num_modes", "1", "--cpu", "6", "--verbosity", "1"])

    jobs = [
        ("lkina_metal", lkina_cmd(["--metal_mode", mode], str(DOCK_DIR / f"{tag}_lkina_{mode}.pdbqt"))),
        ("ad4_std",     lkina_cmd(["--no_auto_metal"], str(DOCK_DIR / f"{tag}_ad4.pdbqt"))),
        ("vina127",     [VINA127, "--receptor", str(rec_qt), "--ligand", str(lig_qt)] + box +
                        ["--out", str(DOCK_DIR / f"{tag}_vina.pdbqt"), "--seed", "42",
                         "--exhaustiveness", "8", "--num_modes", "1", "--cpu", "6"]),
    ]
    for label, cmd in jobs:
        outp = Path(cmd[cmd.index("--out") + 1])
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            ok = outp.exists() and r.returncode == 0
            errtail = (r.stderr or "").strip().splitlines()[-1:] if (r.stderr or "").strip() else []
        except subprocess.TimeoutExpired:
            ok, errtail = False, ["timeout"]
        if ok:
            rmsd = rmsd_best_effort(tag, str(outp))
            try:
                eline = next(ln for ln in open(outp) if "VINA RESULT" in ln)
                energy = float(eline.split()[3])
            except Exception:
                energy = None
            res[label] = {"rmsd": rmsd, "energy": energy}
        else:
            res[label] = {"rmsd": None, "energy": None,
                          "fail": (errtail[0][:120] if errtail else "no output")}
    return res


def main():
    metal_code = sys.argv[1]      # ZN / FE / CU
    want = int(sys.argv[2] if len(sys.argv) > 2 else 40)
    ids, total = rcsb_search(metal_code, 200)
    print(f"[{metal_code}] pool={total}, trying up to {want}")
    done, skipped = [], []
    for pid in ids:
        if len(done) >= want:
            break
        tag_base = f"{pid}_{metal_code}"
        src = PDB_DIR / f"{pid}.pdb"
        if not src.exists():
            r = sh(f'curl -s --max-time 60 -o "{src}" https://files.rcsb.org/download/{pid}.pdb')
            if not src.exists():
                skipped.append(pid); continue
        try:
            built = build_files(pid, metal_code)
        except Exception as ex:
            skipped.append(pid); continue
        if built is None:
            skipped.append(pid)
            # keep the cached pdb (a failed ligand-pick wastes the download);
            # only skip quietly
            continue
        print(f"  [{pid}] prepared, docking...", flush=True)
        try:
            r = run_one(tag_base)
        except subprocess.TimeoutExpired:
            r = {"id": tag_base, "status": "timeout"}
        except Exception as ex:
            r = {"id": tag_base, "status": f"error:{ex.__class__.__name__}"}
        done.append(r)
        print("   ->", json.dumps(r)[:160], flush=True)
    outf = RES_DIR / f"redock_{metal_code.lower()}_{want}.json"
    json.dump({"metal": metal_code, "pool_total": total, "n_attempted": len(done),
               "skipped": len(skipped), "results": done},
              open(outf, "w"), indent=1)
    print("saved:", outf)


if __name__ == "__main__":
    main()
