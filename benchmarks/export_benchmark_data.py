#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export all LKina benchmark test data into standalone, extractable CSV tables.

Each benchmark family becomes one CSV in data_export/ (UTF-8 with BOM so Excel
opens Chinese correctly), plus a README index. Source JSONs are the archived
result files in benchmarks/ (byi/LKina repo).
"""
import json, os, csv, glob, sys

SRC = "/Users/luoxiaowen/Desktop/LKDock/byi/LKina/benchmarks"
OUT = "/Users/luoxiaowen/Desktop/LKDock/LKina论文/data_export"
os.makedirs(OUT, exist_ok=True)

def jload(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)

def write_csv(name, rows, header):
    path = os.path.join(OUT, name)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  {name}: {len(rows)} rows")
    return path

print("=== 1/10 110 metal modes ===")
rows = []
for r in jload(f"{SRC}/metal_coverage_results_all.json"):
    dd = r.get("d_lig")
    d0 = r.get("d0")
    err = round(abs(dd - d0), 3) if (dd is not None and d0) else None
    rows.append({
        "mode": r.get("mode"), "token": r.get("token"), "pseudo": r.get("pseudo"),
        "d0_A": d0, "d_lig_A": round(dd, 3) if dd else None, "abs_err_A": err,
        "e_ideal_kcal": r.get("e_ideal"), "e_far_kcal": r.get("e_far"),
        "well_exists": "YES" if (r.get("e_ideal") is not None and r.get("e_far") is not None and r["e_ideal"] < r["e_far"] + 0.3) else "NO",
        "lkina_rc": r.get("lkina_rc"), "vina_rc": r.get("vina_rc"), "vina_d_A": r.get("vina_d"),
    })
write_csv("metal_modes_110.csv", rows, [
    "mode","token","pseudo","d0_A","d_lig_A","abs_err_A","e_ideal_kcal","e_far_kcal",
    "well_exists","lkina_rc","vina_rc","vina_d_A"])

print("=== 2/10 4JC engine comparison ===")
rows = []
for r in jload(f"{SRC}/real_systems/4JC_engine_comparison.json"):
    rows.append({
        "engine": r.get("engine"), "energy_kcal": r.get("energy"),
        "d_zn_A": r.get("d_zn"), "d_zn_NA_A": r.get("d_zn_NA"),
        "nearest_atom": r.get("nearest"), "coord_report": r.get("coord"),
        "energy_source": r.get("energy_source"),
    })
write_csv("4JC_comparison.csv", rows, [
    "engine","energy_kcal","d_zn_A","d_zn_NA_A","nearest_atom","coord_report","energy_source"])

print("=== 3/10 metallocomplex full pool (n=20) ===")
rows = []
for mc in ["pt", "pd", "ru", "os", "re"]:
    d = jload(f"{SRC}/metallocomplex_results/metallocomplex_{mc}.json")
    for r in d.get("results", []):
        lk = r.get("lkina_metal", {}) or {}
        ad = r.get("ad4_std", {}) or {}
        vn = r.get("vina127", {}) or {}
        rows.append({
            "pdb": r.get("id"), "metal": r.get("metal"), "ligand_resn": r.get("lig"),
            "lkina_metal_rmsd_A": lk.get("rmsd"), "lkina_metal_energy": lk.get("energy"),
            "lkina_metal_geom_penalty": lk.get("metal_geom"),
            "ad4_std_rmsd_A": ad.get("rmsd"), "ad4_std_energy": ad.get("energy"),
            "vina127_rmsd_A": vn.get("rmsd"), "vina127_fail": (vn.get("fail") or "")[:60],
        })
write_csv("metallocomplex_pool.csv", rows, [
    "pdb","metal","ligand_resn","lkina_metal_rmsd_A","lkina_metal_energy",
    "lkina_metal_geom_penalty","ad4_std_rmsd_A","ad4_std_energy","vina127_rmsd_A","vina127_fail"])

print("=== 4/10 C3 deep-pocket ablation ===")
rows = []
d = jload(f"{SRC}/c3_ablation_results.json")
for preset, info in d.get("results", {}).items():
    for variant, v in info.get("variants", {}).items():
        rows.append({
            "preset": preset, "r0_A": info.get("r0"), "theta_deg": info.get("theta"),
            "variant": variant, "rc": v.get("rc"),
            "reactive_dist_A": v.get("REACTIVE_DIST"), "reactive_angle_deg": v.get("REACTIVE_ANGLE"),
            "nac": v.get("REACTIVE_NAC"), "energy_kcal": v.get("energy"),
            "converged": v.get("converged"),
        })
write_csv("c3_ablation.csv", rows, [
    "preset","r0_A","theta_deg","variant","rc","reactive_dist_A",
    "reactive_angle_deg","nac","energy_kcal","converged"])

print("=== 5/10 metal_soft_weight ablation ===")
rows = []
d = jload(f"{SRC}/soft_weight_results.json")
for sysname, info in d.items():
    for w, v in info.get("w", {}).items():
        rows.append({
            "system": sysname, "metal_mode": info.get("mode"), "w": w,
            "rc": v.get("rc"), "donor_metal_A": v.get("donor_metal"),
            "metal_geo_e_kcal": v.get("metal_geo_e"), "energy_kcal": v.get("energy"),
            "fail": v.get("fail"),
        })
write_csv("soft_weight_ablation.csv", rows, [
    "system","metal_mode","w","rc","donor_metal_A","metal_geo_e_kcal","energy_kcal","fail"])

print("=== 6/10 covalent presets (6 presets) ===")
rows = []
d = jload(f"{SRC}/covalent_full_results.json")
for p in d.get("presets", []):
    p4 = p.get("p4_vdw_scale", {}) or {}
    c3 = p.get("c3_two_step", {}) or {}
    p1 = p.get("p1", {}) or {}
    p12 = p.get("p12", {}) or {}
    rows.append({
        "preset": p.get("preset"), "r0_A": p.get("r0"), "theta_deg": p.get("theta"),
        "p1_rc": p1.get("rc"), "p1_nac": p1.get("nac"), "p1_dist_A": p1.get("dist"),
        "p12_rc": p12.get("rc"), "p12_nac": p12.get("nac"),
        "p12_dist_A": p12.get("dist"), "p12_angle_deg": p12.get("angle"),
        "p4_scale0_energy": (p4.get("0.0") or {}).get("energy"),
        "p4_scale05_energy": (p4.get("0.5") or {}).get("energy"),
        "p4_scale10_energy": (p4.get("1.0") or {}).get("energy"),
        "c3_rc": c3.get("rc"), "c3_pose_count": c3.get("size"),
    })
write_csv("covalent_presets.csv", rows, [
    "preset","r0_A","theta_deg","p1_rc","p1_nac","p1_dist_A","p12_rc","p12_nac",
    "p12_dist_A","p12_angle_deg","p4_scale0_energy","p4_scale05_energy",
    "p4_scale10_energy","c3_rc","c3_pose_count"])

print("=== 7/10 feature family (pseudoatom / BVS / water bridge / metal-as-ligand) ===")
rows = []
d = jload(f"{SRC}/feature_family_results.json")
for r in d.get("pseudoatom_geometry", []):
    rows.append({
        "family": "pseudoatom_geometry", "mode": r.get("mode"), "pseudo": r.get("pseudo"),
        "d0_A": r.get("d0"), "donors_checked": r.get("n_donors_checked"),
        "ideal_count": r.get("n_ideal"), "vacancies": r.get("vacancies"), "rc": r.get("rc"),
        "note": " ".join(map(str, r.get("dists", [])))[:80],
    })
for r in d.get("bvs_inference", []):
    rows.append({
        "family": "bvs_inference", "mode": r.get("expected"),
        "expected": r.get("expected"), "detected": r.get("detected"),
        "correct": "YES" if r.get("ok") else "NO",
        "note": r.get("label"),
    })
wb = d.get("water_bridge", {})
if isinstance(wb, dict):
    rows.append({
        "family": "water_bridge", "mode": wb.get("config"), "rc": wb.get("rc"),
        "note": f"metal_water_e={wb.get('metal_water_e')} metal_geo_e={wb.get('metal_geo_e')} coord={wb.get('metal_coord')}",
    })
for r in d.get("metal_as_ligand", []):
    rows.append({
        "family": "metal_as_ligand", "mode": r.get("label"), "rc": r.get("rc"),
        "note": f"angle={r.get('angle')} geom={r.get('ligand_metal_geom')} site={r.get('ligand_metal_site','')[:60]}",
    })
write_csv("feature_family.csv", rows, [
    "family","mode","pseudo","d0_A","donors_checked","ideal_count","vacancies",
    "expected","detected","correct","rc","note"])

print("=== 8/10 reactive presets (P1 vs P1+P2) ===")
rows = []
for r in jload(f"{SRC}/reactive_presets_results.json"):
    rows.append({
        "preset": r.get("preset"), "r0_A": r.get("r0"), "theta_deg": r.get("theta"),
        "p1_rc": r.get("p1_rc"), "p1_nac": r.get("p1_nac"), "p1_dist_A": r.get("p1_dist"),
        "p12_rc": r.get("p12_rc"), "p12_nac": r.get("p12_nac"),
        "p12_dist_A": r.get("p12_dist"), "p12_angle_deg": r.get("p12_angle"),
        "p12_energy_kcal": r.get("p12_re"),
    })
write_csv("reactive_presets.csv", rows, [
    "preset","r0_A","theta_deg","p1_rc","p1_nac","p1_dist_A",
    "p12_rc","p12_nac","p12_dist_A","p12_angle_deg","p12_energy_kcal"])

print("=== 9/10 redock summary (per metal, per engine) ===")
rows = []
d = jload(f"{SRC}/redock_benchmark/results/redock_summary.json")
for metal, m in d.items():
    for engine in ["lkina_metal", "ad4_std", "vina127"]:
        e = m.get(engine, {}) or {}
        rows.append({
            "metal": metal, "n_attempted": m.get("n_attempted"), "n_docked": m.get("n_docked"),
            "engine": engine, "n_rmsd": e.get("n_rmsd"), "success_le2A": e.get("success"),
            "success_pct": e.get("rate_pct"), "mean_rmsd_A": e.get("mean_rmsd"),
        })
write_csv("redock_summary.csv", rows, [
    "metal","n_attempted","n_docked","engine","n_rmsd","success_le2A","success_pct","mean_rmsd_A"])

print("=== 10/10 redock donor-metal distance summary ===")
rows = []
d = jload(f"{SRC}/redock_benchmark/results/donor_metal_distance_summary.json")
for engine, e in d.items():
    rows.append({
        "engine": engine, "n": e.get("n"), "mean_A": e.get("mean"),
        "median_A": e.get("median"), "le_3A_count": e.get("le_3A"),
        "le_3A_pct": round(100.0 * e.get("le_3A", 0) / e.get("n", 1), 1) if e.get("n") else None,
    })
write_csv("redock_donor_metal.csv", rows, [
    "engine","n","mean_A","median_A","le_3A_count","le_3A_pct"])

print("=== 汇总 metallocomplex_summary.json (full pool) ===")
summ = {"n_total": 0, "by_metal": {}}
for mc in ["PT", "PD", "RU", "OS", "RE"]:
    d = jload(f"{SRC}/metallocomplex_results/metallocomplex_{mc.lower()}.json")
    ok = [r for r in d.get("results", []) if r.get("status") == "ok"]
    b = {"n": len(ok), "lk_rmsd_ok": 0, "ad_rmsd_ok": 0, "vn_rmsd_ok": 0, "vn_fail": 0, "lk_geom": []}
    for r in ok:
        lk, ad, vn = r.get("lkina_metal", {}), r.get("ad4_std", {}), r.get("vina127", {})
        if lk.get("rmsd") is not None and lk["rmsd"] <= 2.0: b["lk_rmsd_ok"] += 1
        if ad.get("rmsd") is not None and ad["rmsd"] <= 2.0: b["ad_rmsd_ok"] += 1
        if vn.get("rmsd") is not None and vn["rmsd"] <= 2.0: b["vn_rmsd_ok"] += 1
        if vn.get("fail"): b["vn_fail"] += 1
        if lk.get("metal_geom") is not None: b["lk_geom"].append(lk["metal_geom"])
    summ["by_metal"][mc] = b
    summ["n_total"] += len(ok)
with open(f"{SRC}/metallocomplex_results/metallocomplex_summary.json", "w", encoding="utf-8") as f:
    json.dump(summ, f, ensure_ascii=False, indent=1)
print(f"  metallocomplex_summary.json updated: n_total={summ['n_total']}")

print("\nAll CSV exported to", OUT)
