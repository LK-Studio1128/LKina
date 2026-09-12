#!/usr/bin/env python3
"""Aggregate /tmp/topk_results.jsonl -> v1.0.2 redock benchmark tables.

Schema-compatible with the archived files:
  - results/redock_{zn,fe,cu}_120.json (+ root copies)
        per-system per-engine {rmsd (top-1), energy, rmsds [9], dmin, n_models}
  - results/redock_summary.json        per-metal {engine: rate_pct...} + top3/top5
  - results/donor_metal_distance_summary.json
  - results/topk_seed_summary.json     multiseed stats (5 seeds)
  - /tmp/topk_per_system.csv
"""
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

BASE = Path("/Users/luoxiaowen/Desktop/LKDock/byi/LKina/benchmarks/redock_benchmark")
MANU = Path("/Users/luoxiaowen/Desktop/LKDock/LKina论文/benchmarks/redock_benchmark")
RES_DIR = MANU / "results"
JSONL = Path("/tmp/topk_results.jsonl")
DOCK2 = BASE / "docking_v102"

DONOR_TYPES = {"NA", "NS", "N", "OA", "OS", "O", "SA", "S"}
ENGINES = ["lkina_metal", "ad4_std", "vina127"]
SEEDS = [42, 7, 2026, 3407, 12345]


def donor_min_dist(pose_text, metals):
    dmin = None
    for ln in pose_text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        adt = ln[77:79].strip().upper()
        if adt not in DONOR_TYPES:
            continue
        x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        for mx, my, mz in metals:
            d = math.dist((x, y, z), (mx, my, mz))
            if dmin is None or d < dmin:
                dmin = d
    return dmin


def first_model_text(txt):
    parts = txt.split("MODEL")
    return parts[1].split("ENDMDL", 1)[0] if len(parts) > 1 else txt


def main():
    runs = [json.loads(l) for l in open(JSONL)]
    idx = {(r["tag"], r["engine"], r["seed"]): r for r in runs if r["status"] == "ok"}
    engines = ENGINES
    metals_of = {}
    for m in ["zn", "fe", "cu"]:
        d = json.load(open(BASE / f"redock_{m}_120.json"))
        for r in d["results"]:
            if r.get("status") == "ok":
                metals_of[r["id"]] = m.upper()
    tags = sorted(metals_of)
    n42 = sum(1 for t in tags for e in engines if (t, e, 42) in idx)
    print(f"runs ok: {len(idx)} | systems: {len(tags)} | seed-42 engine-jobs: {n42}/312")
    if n42 < 312:
        print("WARNING: seed-42 batch incomplete; aggregating partial data")

    out120 = {m: [] for m in ["ZN", "FE", "CU"]}
    dmins = {e: {} for e in engines}
    csv = ["tag,metal,engine,top1_rmsd,top2_rmsd,top3_rmsd,top5_rmsd,dmin_best_A,best_E,n_models"]
    for tag in tags:
        metal = metals_of[tag]
        meta = json.load(open(RES_DIR / f"{tag}_meta.json"))
        metals = [tuple(m) for m in meta["metals"]]
        entry = {"id": tag, "status": "ok", "lig": meta.get("lig_name", "")}
        for e in engines:
            r = idx.get((tag, e, 42))
            if not r:
                entry[e] = {"rmsd": None, "energy": None, "fail": "v1.0.2 rerun: not completed"}
                continue
            outp = DOCK2 / f"{tag}_{e}_s42.pdbqt"
            dmin = None
            if outp.exists():
                first = first_model_text(open(outp, errors="ignore").read())
                dmin = donor_min_dist(first, metals)
                dmin = round(dmin, 2) if dmin is not None else None
            rs = r["rmsds"]
            entry[e] = {"rmsd": rs[0] if rs else None,
                        "energy": r["energies"][0] if r["energies"] else None,
                        "rmsds": rs, "dmin": dmin, "n_models": r["n_models"]}
            if dmin is not None:
                dmins[e][tag] = dmin
            csv.append(f"{tag},{metal},{e},{rs[0] if rs else ''},"
                       f"{rs[1] if len(rs) > 1 else ''},{rs[2] if len(rs) > 2 else ''},"
                       f"{rs[4] if len(rs) > 4 else ''},{dmin if dmin is not None else ''},"
                       f"{r['energies'][0] if r['energies'] else ''},{r['n_models']}")
        out120[metal].append(entry)

    for m in ["ZN", "FE", "CU"]:
        old = json.load(open(BASE / f"redock_{m.lower()}_120.json"))
        old["results"] = out120[m]
        old["source"] = ("v1.0.2 rerun 2026-09-12 (includes 0b164d2 TZ-well fix; "
                         "exh 8, cpu 6, num_modes 9, seed 42 shown; 5-seed set archived)")
        for dest in [BASE / f"redock_{m.lower()}_120.json",
                     BASE / "results" / f"redock_{m.lower()}_120.json",
                     RES_DIR / f"redock_{m.lower()}_120.json"]:
            json.dump(old, open(dest, "w"), ensure_ascii=False, indent=1)

    # ---------- per-metal summary (compat + top-k) ----------
    # Semantics: r["rmsds"] is in Vina OUTPUT order (rank 1 = best-energy pose).
    #   top-1 = rmsds[0];  top-k = min(rmsds[:k]) over the first k ranked poses.
    def rate(engine, metal, thr, k):
        vals = []
        for tag in tags:
            if metals_of[tag] != metal:
                continue
            r = idx.get((tag, engine, 42))
            if not r or not r["rmsds"]:
                continue
            vals.append(min(r["rmsds"][:k]) <= thr)
        n = len(vals)
        return (sum(vals), n)

    redock_summary = {}
    for metal in ["ZN", "FE", "CU"]:
        attempted = sum(1 for t in tags if metals_of[t] == metal)
        redock_summary[metal] = {"n_attempted": attempted}
        for e in engines:
            s1, n = rate(e, metal, 2.0, 1)
            s3, _ = rate(e, metal, 2.0, 3)
            s5, _ = rate(e, metal, 2.0, 5)
            g1, gn = rate(e, metal, 3.0, 1)
            rms = [idx[(t, e, 42)]["rmsds"][0] for t in tags
                   if metals_of[t] == metal and (t, e, 42) in idx
                   and idx[(t, e, 42)]["rmsds"]]
            redock_summary[metal][e] = {
                "n_rmsd": n, "success": s1, "rate_pct": round(100 * s1 / n, 1) if n else 0.0,
                "mean_rmsd": round(statistics.mean(rms), 2) if rms else None,
                "top3_le2": s3, "top5_le2": s5, "le3_top1": g1,
                "rate_top3_pct": round(100 * s3 / n, 1) if n else 0.0,
                "rate_top5_pct": round(100 * s5 / n, 1) if n else 0.0,
            }

    dsum = {}
    for e in engines:
        vals = list(dmins[e].values())
        dsum[e] = {"n": len(vals),
                   "mean": round(statistics.mean(vals), 2),
                   "median": round(statistics.median(vals), 2),
                   "le_3A": sum(v <= 3.0 for v in vals)} if vals else {"n": 0}

    # ---------- multiseed stats (5 seeds) ----------
    multi = {}
    for e in engines:
        rows = []
        for tag in tags:
            es, rs = [], []
            for s in SEEDS:
                r = idx.get((tag, e, s))
                if r and r["energies"]:
                    es.append(r["energies"][0])
                if r and r["rmsds"]:
                    rs.append(r["rmsds"][0])
            if len(es) == 5 and len(rs) == 5:
                rows.append({"tag": tag, "metal": metals_of[tag],
                             "E_mean": round(statistics.mean(es), 2),
                             "E_std": round(statistics.stdev(es), 3),
                             "rmsd_top1_mean": round(statistics.mean(rs), 2),
                             "rmsd_top1_std": round(statistics.stdev(rs), 2),
                             "success_all5": all(x <= 2.0 for x in rs),
                             "success_any5": any(x <= 2.0 for x in rs),
                             "le3_all5": None, "dmin_std": None})
        # donor<=3A stability across seeds for subset systems
        per_seed_le3 = {s: 0 for s in SEEDS}
        per_seed_tot = 0
        dmin_by_tag_seed = {}
        for tag in tags:
            for s in SEEDS:
                r = idx.get((tag, e, s))
                if not r:
                    continue
                outp = DOCK2 / f"{tag}_{e}_s{s}.pdbqt"
                if not outp.exists():
                    continue
                meta = json.load(open(RES_DIR / f"{tag}_meta.json"))
                dmin = donor_min_dist(first_model_text(open(outp, errors="ignore").read()),
                                      [tuple(m) for m in meta["metals"]])
                if dmin is None:
                    continue
                dmin_by_tag_seed[(tag, s)] = dmin
        per_seed_tot = len({t for (t, s) in dmin_by_tag_seed})
        for s in SEEDS:
            per_seed_le3[s] = sum(1 for (t, ss), d in dmin_by_tag_seed.items()
                                  if ss == s and d <= 3.0)
        e_tot_sys = per_seed_tot
        multi[e] = {"n_systems_5seed": len(rows),
                    "E_std_median": round(statistics.median([r["E_std"] for r in rows]), 3) if rows else None,
                    "rmsd_top1_std_median": round(statistics.median([r["rmsd_top1_std"] for r in rows]), 2) if rows else None,
                    "success_consistent_all5": sum(r["success_all5"] for r in rows),
                    "success_ever5": sum(r["success_any5"] for r in rows),
                    "donor_le3_per_seed": per_seed_le3,
                    "donor_le3_systems": e_tot_sys}
        json.dump(rows, open(BASE / "results" / f"topk_persystem_{e}.json", "w"),
                  ensure_ascii=False, indent=1)

    summary = {"protocol": "exh 8, cpu 6, num_modes 9, box = ligand span + 12 A (min 10), Open Babel prep",
               "binary": "LKina v1.0.2 (build/mac/release, includes 0b164d2 TZ-well fix); Vina 1.2.7 (PLIP build)",
               "date": "2026-09-12", "seeds": SEEDS,
               "per_metal": redock_summary, "donor_summary_seed42": dsum,
               "multiseed": multi}
    # write everywhere the figure/export scripts and the repo expect them
    for dest in [BASE / "results" / "topk_seed_summary.json", BASE / "topk_seed_summary.json"]:
        json.dump(summary, open(dest, "w"), ensure_ascii=False, indent=1)
    for dest in [RES_DIR / "donor_metal_distance_summary.json",
                 BASE / "results" / "donor_metal_distance_summary.json",
                 BASE / "donor_metal_distance_summary.json"]:
        json.dump(dsum, open(dest, "w"), indent=1)
    for dest in [RES_DIR / "redock_summary.json",
                 BASE / "results" / "redock_summary.json",
                 BASE / "redock_summary.json"]:
        json.dump(redock_summary, open(dest, "w"), indent=1)
    Path("/tmp/topk_per_system.csv").write_text("\n".join(csv) + "\n")

    print("== per-metal top-1 RMSD<=2.0 (seed42) ==")
    for metal in ["ZN", "FE", "CU"]:
        row = " | ".join(f"{e}: {redock_summary[metal][e]['rate_pct']}% (top3 {redock_summary[metal][e]['rate_top3_pct']}%)" for e in engines)
        print(f" {metal}: {row}")
    print("== donor<=3A (seed42, engine-dockable subset) ==")
    for e in engines:
        print(f" {e}: {dsum[e]['le_3A']}/{dsum[e]['n']} ({round(100*dsum[e]['le_3A']/max(1,dsum[e]['n']),1)}%) median {dsum[e]['median']}")
    print("written.")


if __name__ == "__main__":
    main()
