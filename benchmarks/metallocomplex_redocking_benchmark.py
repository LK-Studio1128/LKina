#!/usr/bin/env python3
"""
LKina metallocomplex redocking benchmark — "metal as ligand" (P1).

Protocol (MetalDock-comparable, per paper Section 2.8):
  1. Download PDB entries where a coordination-metal HET ligand
     (PT / PD / RU / OS / RE) is the bound ligand (metal complex as ligand).
  2. Extract receptor (protein, no co-crystallized ligand) and the
     metal-complex ligand (HETATM residue containing the metal).
  3. Prepare PDBQT: receptor via obabel rigid; ligand via obabel.
  4. Redock with:
       a) LKina AD4 (auto-detected metal-as-ligand reverse pair potentials
          + `--ligand_metal_geometry_weight 1.0` ligand-side geometry QC)
       b) LKina AD4 standard (`--no_auto_metal`)  -- no metal-aware terms
       c) Vina 1.2.7 (baseline; expected to fail on metal atom types)
  5. Metrics per system:
       - top-1 RMSD <= 2.0 A vs crystal pose (redock success)
       - LIGAND_METAL_GEOM penalty from the docked pose (geometry QC active)
  6. Aggregate by metal family (Pt/Pd/Ru/Os/Re).

Answers the MetalDock critique: do LKina's reverse pair potentials help
redock metal-complex ligands back into their protein pockets?

Usage:
  python3 benchmarks/metallocomplex_redocking_benchmark.py PT [N]
  python3 benchmarks/metallocomplex_redocking_benchmark.py --all
"""
import json, math, os, subprocess, sys, time, urllib.request
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).parent
PDB_DIR = BASE / "metallocomplex_pdb"
PREP_DIR = BASE / "metallocomplex_prepared"
DOCK_DIR = BASE / "metallocomplex_docking"
RES_DIR = BASE / "metallocomplex_results"
for d in (PDB_DIR, PREP_DIR, DOCK_DIR, RES_DIR):
    d.mkdir(exist_ok=True)

LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
VINA127 = "/Volumes/LK/软件资料/软件开发/分子对接全流程工具-PLIP版/vina"

METAL_CODES = ["PT", "PD", "RU", "OS", "RE"]

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
         "HSD HSE HSP MSE SEP TPO PTR CSO CME OCS CAS MHO SEC PYL".split())
SKIP_RES = {"HOH", "DOD", "WAT", "NA", "CL", "K", "MG", "CA", "SO4", "PO4", "GOL",
            "EDO", "ACT", "NO3", "ZN", "FE", "CU", "MN", "CO", "NI", "CD", "HG"}


def sh(cmd, cwd=None, timeout=600):
    return subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True, timeout=timeout)


def rcsb_search(metal_code, max_results=100):
    """Entries with a non-polymer HET residue containing the metal (as ligand)."""
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                            "operator": "exact_match", "value": metal_code}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "less_or_equal", "value": 2.8}}]},
        "request_options": {"paginate": {"start": 0, "rows": min(max_results, 500)},
                            "results_verbosity": "compact"},
        "return_type": "entry"}
    req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",
                                 data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
        return (d.get("identifier_set") or d.get("result_set") or []), int(d.get("total_count") or 0)
    except Exception as e:
        print("  RCSB search failed:", e)
        return [], 0


def parse_pdb(path):
    atoms = []
    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                resn = line[17:20].strip()
                chain = line[21:22].strip()
                resi = line[22:26].strip()
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                elem = line[76:78].strip().upper()
                if not elem:
                    elem = "".join(c for c in name if c.isalpha())[:1].upper()
                atoms.append({"rec": line[:6].strip(), "name": name, "resn": resn,
                              "chain": chain, "resi": resi, "xyz": xyz, "elem": elem, "line": line})
    return atoms


def dist(a, b):
    return math.dist(a["xyz"], b["xyz"])


def select_metal_ligand(atoms, metal_code):
    """Largest HETATM residue containing the metal atom, excluding water/ions."""
    metals = [a for a in atoms if a["elem"] == metal_code and a["rec"] == "HETATM"]
    if not metals:
        return None
    het_res = defaultdict(list)
    for a in atoms:
        if a["resn"] not in AA and a["resn"] not in SKIP_RES and a["rec"] == "HETATM":
            het_res[(a["chain"], a["resi"], a["resn"])].append(a)
    best, best_n = None, 0
    for key, group in het_res.items():
        if not any(a["elem"] == metal_code for a in group):
            continue
        n_heavy = sum(1 for a in group if a["elem"] != "H")
        if n_heavy < 8:  # metal complex must be a real ligand, not a bare ion
            continue
        if n_heavy > best_n:
            best, best_n = group, n_heavy
    return best


def build_files(pdb_id, metal_code):
    src = PDB_DIR / f"{pdb_id}.pdb"
    atoms = parse_pdb(src)
    lig = select_metal_ligand(atoms, metal_code)
    if lig is None:
        return None
    lig_ids = {(a["chain"], a["resi"], a["resn"]) for a in lig}
    rec_lines, seen = [], set()
    for a in atoms:
        lid = (a["chain"], a["resi"], a["resn"])
        if lid in lig_ids:
            continue
        if a["resn"] in SKIP_RES or a["resn"] in ("HOH", "DOD", "WAT"):
            continue
        if a["resn"] in AA or a["rec"] == "ATOM":
            alt = a["line"][16]
            key = (a["chain"], a["resi"], a["name"], a["resn"])
            if alt not in (" ", "A") or key in seen:
                continue
            seen.add(key)
            rec_lines.append(a["line"])
    out_lig, seenl = [], set()
    for a in lig:
        alt = a["line"][16]
        if alt not in (" ", "A") or a["name"] in seenl:
            continue
        seenl.add(a["name"])
        out_lig.append(a["line"])

    tag = f"{pdb_id}_{metal_code}"
    rp = PREP_DIR / f"{tag}_rec.pdb"
    lp = PREP_DIR / f"{tag}_lig.pdb"
    (PREP_DIR).mkdir(exist_ok=True)
    rp.write_text("\n".join(ln.rstrip() for ln in rec_lines if ln.strip()) + "\nEND\n")
    lp.write_text("\n".join(ln.rstrip() for ln in out_lig if ln.strip()) + "\nTER\nEND\n")

    xyz = [(float(ln[30:38]), float(ln[38:46]), float(ln[46:54])) for ln in out_lig]
    meta = {"pdb": pdb_id, "metal": metal_code,
            "lig_name": lig[0]["resn"], "n_heavy": sum(1 for a in lig if a["elem"] != "H"),
            "lig_center": [round(sum(p[i] for p in xyz) / len(xyz), 3) for i in range(3)],
            "span": [round(max(p[i] for p in xyz) - min(p[i] for p in xyz), 1) for i in range(3)]}
    (RES_DIR / f"{tag}_meta.json").write_text(json.dumps(meta))
    return lp, rp, xyz


def prep_pdbqt(tag):
    lp = PREP_DIR / f"{tag}_lig.pdb"
    rp = PREP_DIR / f"{tag}_rec.pdb"
    lig_qt = PREP_DIR / f"{tag}_lig.pdbqt"
    rec_qt = PREP_DIR / f"{tag}_rec.pdbqt"
    if not lig_qt.exists() or lig_qt.stat().st_size == 0:
        sh(f'obabel "{lp}" -O "{lig_qt}" 2>/dev/null')
        if not lig_qt.exists() or lig_qt.stat().st_size == 0:
            return None
        txt = lig_qt.read_text()
        if "ROOT" not in txt or txt.count("ROOT") > 1:
            return None
    if not rec_qt.exists() or rec_qt.stat().st_size == 0:
        sh(f'obabel "{rp}" -O "{rec_qt}" -xr 2>/dev/null')
        if not rec_qt.exists() or rec_qt.stat().st_size == 0:
            return None
        lines = [ln for ln in rec_qt.read_text().splitlines()
                 if ln.strip() not in ("ROOT", "ENDROOT", "TORSDOF 0", "BRANCH", "ENDBRANCH")]
        rec_qt.write_text("\n".join(lines) + "\n")
    return rec_qt, lig_qt


_AD_ELEM = {"A": "C", "OA": "O", "NA": "N", "N": "N", "NS": "N", "SA": "S", "S": "S",
            "OS": "O", "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P",
            "SI": "Si", "B": "B", "SE": "Se", "PT": "Pt", "PD": "Pd", "RU": "Ru",
            "OS": "Os", "RE": "Re", "IR": "Ir", "RH": "Rh", "AU": "Au", "AG": "Ag"}


def _heavy_atoms_pdbqt(path):
    out = []
    for ln in open(path, errors="ignore"):
        if ln.startswith(("ATOM", "HETATM")):
            adt = ln[77:79].strip().upper()
            name = ln[12:16].strip()
            elem = _AD_ELEM.get(adt) or "".join(c for c in name if c.isalpha())[:1].upper()
            if elem.startswith("H"):
                continue
            out.append((elem, (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))))
    return out


def rmsd_best_effort(tag, docked_pdbqt):
    try:
        ref = _heavy_atoms_pdbqt(PREP_DIR / f"{tag}_lig.pdbqt")
    except FileNotFoundError:
        return None
    txt = open(docked_pdbqt, errors="ignore").read()
    models = txt.split("MODEL")
    seg = models[1] if len(models) > 1 else txt
    tmp = Path("/tmp") / f"_mal_{tag}.pdbqt"
    tmp.write_text(seg)
    dock = _heavy_atoms_pdbqt(tmp)
    tmp.unlink(missing_ok=True)
    if not ref or not dock or len(ref) != len(dock):
        return None
    if {e for e, _ in ref} != {e for e, _ in dock}:
        return None
    s = sum(math.dist(r, d) ** 2 for (_, r), (_, d) in zip(ref, dock)) / len(ref)
    return round(math.sqrt(s), 2)


def lig_metal_penalty(docked_pdbqt):
    """Parse LIGAND_METAL_GEOM / LIGAND_METAL_SITE from the docked pose."""
    geom, site = None, None
    for ln in open(docked_pdbqt, errors="ignore"):
        if "LIGAND_METAL_GEOM" in ln:
            try:
                geom = float(ln.split(":")[-1].strip())
            except Exception:
                pass
        if "LIGAND_METAL_SITE" in ln:
            site = ln.strip()[:80]
    return geom, site


def run_one(tag):
    prepped = prep_pdbqt(tag)
    if prepped is None:
        return {"id": tag, "status": "prep_fail"}
    rec_qt, lig_qt = prepped
    meta = json.load(open(RES_DIR / f"{tag}_meta.json"))
    cx, cy, cz = meta["lig_center"]
    sx, sy, sz = [max(s + 8.0, 14) for s in meta["span"]]
    box = ["--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", f"{sx:.0f}", "--size_y", f"{sy:.0f}", "--size_z", f"{sz:.0f}"]
    res = {"id": tag, "status": "ok", "lig": meta["lig_name"], "metal": meta["metal"]}

    def lkina_cmd(extra, out):
        return ([LKINA, "--scoring", "ad4", "--generate_maps"] + extra +
                ["--receptor", str(rec_qt), "--ligand", str(lig_qt)] + box +
                ["--out", out, "--seed", "42", "--exhaustiveness", "8",
                 "--num_modes", "1", "--cpu", "6", "--verbosity", "1"])

    jobs = [
        ("lkina_metal",
         lkina_cmd(["--ligand_metal_geometry_weight", "1.0"],
                   str(DOCK_DIR / f"{tag}_lkina_metal.pdbqt"))),
        ("ad4_std",
         lkina_cmd(["--no_auto_metal"], str(DOCK_DIR / f"{tag}_ad4.pdbqt"))),
        ("vina127",
         [VINA127, "--receptor", str(rec_qt), "--ligand", str(lig_qt)] + box +
         ["--out", str(DOCK_DIR / f"{tag}_vina.pdbqt"), "--seed", "42",
          "--exhaustiveness", "8", "--num_modes", "1", "--cpu", "6"]),
    ]
    for label, cmd in jobs:
        outp = Path(cmd[cmd.index("--out") + 1])
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            ok = outp.exists() and outp.stat().st_size > 0 and r.returncode == 0
            errtail = ((r.stderr or "").strip().splitlines()[-1:] if (r.stderr or "").strip() else [])
        except subprocess.TimeoutExpired:
            ok, errtail = False, ["timeout"]
        if ok:
            entry = {"rmsd": rmsd_best_effort(tag, str(outp))}
            try:
                eline = next(ln for ln in open(outp) if "VINA RESULT" in ln)
                entry["energy"] = float(eline.split()[3])
            except Exception:
                entry["energy"] = None
            geom, site = lig_metal_penalty(str(outp))
            entry["metal_geom"] = geom
            entry["metal_site"] = site
            res[label] = entry
        else:
            res[label] = {"rmsd": None, "energy": None,
                          "fail": (errtail[0][:120] if errtail else "no output")}
    return res


def main():
    args = sys.argv[1:]
    if args and args[0] == "--all":
        codes = METAL_CODES
    elif args:
        codes = [args[0].upper()]
    else:
        codes = METAL_CODES
    want = int(args[1]) if len(args) > 1 else 10

    all_results = []
    for mc in codes:
        print(f"=== metal-as-ligand redocking: {mc} ===")
        ids, total = rcsb_search(mc, 200)
        print(f"  pool={total}, trying up to {want}")
        done, skipped = [], []
        for pid in ids:
            if len(done) >= want:
                break
            src = PDB_DIR / f"{pid}.pdb"
            if not src.exists():
                r = sh(f'curl -s --max-time 60 -o "{src}" https://files.rcsb.org/download/{pid}.pdb')
                if not src.exists():
                    skipped.append(pid)
                    continue
            try:
                built = build_files(pid, mc)
            except Exception:
                skipped.append(pid)
                continue
            if built is None:
                skipped.append(pid)
                continue
            tag = f"{pid}_{mc}"
            print(f"  [{pid}] {built[0].name.split('_')[-2] if False else ''}docking...", flush=True)
            try:
                r = run_one(tag)
            except subprocess.TimeoutExpired:
                r = {"id": tag, "status": "timeout"}
            except Exception as ex:
                r = {"id": tag, "status": f"error:{ex.__class__.__name__}"}
            done.append(r)
            print("   ->", json.dumps(r)[:180], flush=True)
        outf = RES_DIR / f"metallocomplex_{mc.lower()}.json"
        json.dump({"metal": mc, "pool_total": total, "n_attempted": len(done),
                   "skipped": len(skipped), "results": done}, open(outf, "w"), indent=1)
        print("saved:", outf)
        all_results += done

    summary = {"n": len(all_results), "by_metal": defaultdict(dict)}
    for r in all_results:
        if r.get("status") != "ok":
            continue
        m = r.get("metal", "?")
        lk, ad, vn = r.get("lkina_metal", {}), r.get("ad4_std", {}), r.get("vina127", {})
        b = summary["by_metal"].setdefault(m, {"n": 0, "lk_rmsd_ok": 0, "ad_rmsd_ok": 0,
                                               "vn_rmsd_ok": 0, "lk_geom": []})
        b["n"] += 1
        if lk.get("rmsd") is not None and lk["rmsd"] <= 2.0:
            b["lk_rmsd_ok"] += 1
        if ad.get("rmsd") is not None and ad["rmsd"] <= 2.0:
            b["ad_rmsd_ok"] += 1
        if vn.get("rmsd") is not None and vn["rmsd"] <= 2.0:
            b["vn_rmsd_ok"] += 1
        if lk.get("metal_geom") is not None:
            b["lk_geom"].append(lk["metal_geom"])
    (RES_DIR / "metallocomplex_summary.json").write_text(json.dumps(summary, indent=1))
    print("\n=== summary ===")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
