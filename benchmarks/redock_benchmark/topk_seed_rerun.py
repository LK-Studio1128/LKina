#!/usr/bin/env python3
"""Top-k and multi-seed redock rerun (advisor comments 4 & 13).

Same protocol as redock_pipeline.py (--exhaustiveness 8, --cpu 6, same boxes),
but --num_modes 9 and 5 seeds. Per-model RMSD uses the same no-superposition
order-paired method as the original pipeline (rmsd_best_effort).

Usage: python3 topk_seed_rerun.py [tag ...]   (default: 24-system subset)
Output: JSONL (one line per run) appended to /tmp/topk_results.jsonl
"""
import json
import math
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = Path("/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks/redock_benchmark")
PREP_DIR = BASE / "prepared"
RES_DIR = BASE / "results"
OUTDIR = Path("/tmp/topk_docking")
OUTDIR.mkdir(exist_ok=True)
JSONL = Path("/tmp/topk_results.jsonl")

LKINA = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/build/mac/release/LKina"
VINA127 = "/Volumes/LK/软件资料/软件开发/分子对接全流程工具-PLIP版/vina"
METAL_MODES = {"ZN": "zn", "FE": "fe3", "CU": "cu2_jt"}
SEEDS = [42, 7, 2026, 3407, 12345]
BOX_MARGIN = 6.0

SUBSET_24 = ("10TV_ZN 10TX_ZN 10TY_ZN 13MQ_ZN 13NE_ZN 13NI_ZN 13SB_ZN 13SM_ZN "
             "13SQ_ZN 13TB_ZN 1A42_ZN 1A4M_ZN 1A8E_FE 1BSZ_FE 1EOB_FE 1HU9_FE "
             "1LRM_FE 1NO3_FE 1A2V_CU 1D6U_CU 1G3E_CU 1IBY_CU 1JES_CU 1L9P_CU").split()


def build_jobs(tags):
    """Priority order: (1) seed-42 x all systems x engines (table-critical),
    (2) 24-system subset x remaining seeds (variance analysis),
    (3) full cohort x remaining seeds (bonus, keeps running)."""
    jobs = [(t, e, 42) for t in tags for e in ("lkina_metal", "ad4_std", "vina127")]
    jobs += [(t, e, s) for t in SUBSET_24 for s in SEEDS if s != 42
             for e in ("lkina_metal", "ad4_std", "vina127")]
    rest = [t for t in tags if t not in SUBSET_24]
    jobs += [(t, e, s) for t in rest for s in SEEDS if s != 42
             for e in ("lkina_metal", "ad4_std", "vina127")]
    return jobs

_AD_ELEM = {"A": "C", "OA": "O", "NA": "N", "N": "N", "NS": "N", "SA": "S", "S": "S",
            "OS": "O", "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P",
            "SI": "Si", "B": "B", "SE": "Se"}


def _heavy_atoms_pdbqt(text):
    out = []
    for ln in text.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            adt = ln[77:79].strip().upper()
            name = ln[12:16].strip()
            elem = _AD_ELEM.get(adt)
            if elem is None:
                elem = "".join(c for c in name if c.isalpha())[:1].upper()
            if elem.startswith("H"):
                continue
            out.append((elem, (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))))
    return out


def models_of(path):
    """Return list of model segments (text) from a multi-model pdbqt."""
    txt = open(path, errors="ignore").read()
    parts = txt.split("MODEL")
    if len(parts) <= 1:
        return [txt]
    segs = []
    for p in parts[1:]:
        seg = p.split("ENDMDL", 1)[0]
        segs.append(seg)
    return segs


def rmsd_models(path, tag):
    """Per-model no-superposition RMSD vs prepared ligand (order-paired)."""
    ref = _heavy_atoms_pdbqt(open(PREP_DIR / f"{tag}_lig.pdbqt", errors="ignore").read())
    out = []
    for seg in models_of(path):
        dock = _heavy_atoms_pdbqt(seg)
        if not ref or not dock or len(ref) != len(dock):
            out.append(None)
            continue
        if {e for e, _ in ref} != {e for e, _ in dock}:
            out.append(None)
            continue
        s = sum(math.dist(r, d) ** 2 for (_, r), (_, d) in zip(ref, dock)) / len(ref)
        out.append(round(math.sqrt(s), 2))
    return out


def energies_models(path):
    en = []
    for seg in models_of(path):
        for ln in seg.splitlines():
            if "VINA RESULT" in ln:
                en.append(float(ln.split()[3]))
                break
    return en


def run_job(tag, engine, seed):
    meta = json.load(open(RES_DIR / f"{tag}_meta.json"))
    cx, cy, cz = meta["lig_center"]
    sx, sy, sz = [max(s + BOX_MARGIN * 2, 10) for s in meta["span"]]
    metal = meta["metal"]
    box = ["--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", f"{sx:.0f}", "--size_y", f"{sy:.0f}", "--size_z", f"{sz:.0f}"]
    rec_qt, lig_qt = PREP_DIR / f"{tag}_rec.pdbqt", PREP_DIR / f"{tag}_lig.pdbqt"
    out = OUTDIR / f"{tag}_{engine}_s{seed}.pdbqt"
    if engine == "lkina_metal":
        cmd = [LKINA, "--scoring", "ad4", "--generate_maps", "--metal_mode", METAL_MODES[metal],
               "--receptor", str(rec_qt), "--ligand", str(lig_qt)] + box + \
              ["--out", str(out), "--seed", str(seed), "--exhaustiveness", "8",
               "--num_modes", "9", "--cpu", "6", "--verbosity", "1"]
    elif engine == "ad4_std":
        cmd = [LKINA, "--scoring", "ad4", "--generate_maps", "--no_auto_metal",
               "--receptor", str(rec_qt), "--ligand", str(lig_qt)] + box + \
              ["--out", str(out), "--seed", str(seed), "--exhaustiveness", "8",
               "--num_modes", "9", "--cpu", "6", "--verbosity", "1"]
    else:
        cmd = [VINA127, "--receptor", str(rec_qt), "--ligand", str(lig_qt)] + box + \
              ["--out", str(out), "--seed", str(seed), "--exhaustiveness", "8",
               "--num_modes", "9", "--cpu", "6"]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        ok = out.exists() and r.returncode == 0
    except subprocess.TimeoutExpired:
        ok = False
    dt = round(time.time() - t0, 1)
    if not ok:
        return {"tag": tag, "engine": engine, "seed": seed, "status": "fail", "wall_s": dt}
    rmsds = rmsd_models(str(out), tag)
    en = energies_models(str(out))
    return {"tag": tag, "engine": engine, "seed": seed, "status": "ok", "wall_s": dt,
            "rmsds": rmsds, "energies": en, "n_models": len(models_of(str(out)))}


def main():
    if len(sys.argv) > 1:
        tags = sys.argv[1:]
    else:
        allf = Path("/tmp/all_tags.txt")
        if allf.exists():
            tags = open(allf).read().strip().split("\n")[-1].split()
        else:
            tags = SUBSET_24
    jobs = build_jobs(tags)
    # skip jobs already recorded in the JSONL (resume support)
    done = set()
    if JSONL.exists():
        for ln in open(JSONL):
            try:
                r = json.loads(ln)
                done.add((r["tag"], r["engine"], r["seed"]))
            except Exception:
                pass
    todo = [j for j in jobs if j not in done]
    print(f"{len(jobs)} jobs total, {len(jobs)-len(todo)} already done, {len(todo)} to run, 4 workers", flush=True)
    n = len(todo)
    done_n = 0
    with JSONL.open("a") as sink, ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(run_job, *j): j for j in todo}
        for f in as_completed(futs):
            res = f.result()
            sink.write(json.dumps(res) + "\n")
            sink.flush()
            done_n += 1
            if done_n % 12 == 0 or done_n == n:
                print(f"[{done_n}/{n}] last: {res['tag']} {res['engine']} s{res['seed']} "
                      f"{res['status']} {res['wall_s']}s", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
