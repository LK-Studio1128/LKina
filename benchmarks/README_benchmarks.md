# LKina Benchmark Results — Reproducibility Summary

All runs executed on this machine (macOS, Apple Silicon, clang++-built binaries).
Fixed seed 42 throughout. Scripts in this directory are the source of truth.

**Update (round 4, post-engine-fix v1.0.0):** full-corpus re-run with the
±1000 kcal/mol grid clamp and `--no_auto_metal` flag. The live binaries are the
post-fix build; all numbers below come from the v3/v4 scripts.

## Shared PDBQT writer

`pdbqt_util.py` — single source of truth for synthetic receptor/ligand writing.
Absolute-column layout (0-indexed): chain @21, resi @22–25, x/y/z @30–37/38–45/46–53,
occ @54–59, b @60–65, charge @66–76, type @77–78. Verified to parse identically
under LKina and Vina 1.2.7 (earlier column-offset issues resolved).

## 1. Metal-coverage benchmark (110 metal modes — the complete `metal_mode` enum)

Script: `run_all_metal_modes.py` -> `metal_coverage_results_all.json`
Systems: per metal_mode, a synthetic receptor (metal at origin + (n−1)
coordination donors at the ideal distance along the mode geometry) and an
NA-donor probe ligand starting at d0+4 Å.

| Metric | LKina | Vina 1.2.7 |
|---|---|---|
| Docking completed | **110/110** | 34/110 |
| Failed (unsupported atom type, parse time) | 0 | **76** |
| Mean \|d(M–donor) − d0\| (Å) | **0.20** | 0.62 (34 dockable) |
| \|err\| < 0.5 Å | **108/110** | — |
| \|err\| < 1.0 Å | 109/110 | — |
| Coordination well confirmed (E_ideal < E_far) | 104/110 | — |

Per-pseudoatom-class mean |d−d0|: TZ 0.25 (n=2), SQ 0.24 (n=7), MH 0.20 (n=95),
JT 0.19 (n=2), LIN 0.14 (n=4).

Vina 1.2.7 failure mode: `PDBQT parsing error: Atom type <M> is not a valid
AutoDock type (atom types are case-sensitive).` — the 32-type XS system cannot
map LKina's 113-type AD4 atom space; 76/110 metal modes are unparseable.

Vina 1.2.7 success cases (34 runs over: At Ca Cd Co×4 Cu×4 Fe×4 Hg K Mg×2
Mn×5 Na Ni×3 Se U×2 W Zn).

## 1b. Feature-family measurements (`run_feature_family_tests.py` -> `feature_family_results.json`)

- **Pseudoatom geometry check** (`--metal_geometry_check`, five metals):
  zn→TZ 3/3, pt→SQ 3/3, fe→MH 5/5, cu2_jt→JT 4/4, mn3_jt→JT 4/4 donors
  classified ideal/good against the mode-specific nbp r_eq.
- **BVS oxidation-state inference**: **14/14** designed systems correct
  (Fe²⁺/Fe³⁺, Cu⁺/Cu²⁺, Mn²⁺/Mn³⁺, Co²⁺/Co³⁺, V⁴⁺/V⁵⁺, Mo⁴⁺/Mo⁶⁺, Ni²⁺/Ni³⁺),
  auto-detected without `--metal_mode`.
- **Water bridge** (Mg, octahedral 4/6 occupied → 2 vacant sites):
  `METAL_WATER_E −0.000`, `METAL_GEO_E −0.090`,
  `METAL_COORD: Mg donor=NA d=2.28 A`, n_vacant_sites = 2.
- **Metal-as-ligand geometry QC** (Pt(II) square-planar):
  ideal 90° → penalty 0.003 kcal/mol (quality=good);
  distorted 60° variant → penalty 0.903 kcal/mol (≈300×).

## 1c. Covalent framework, four tiers (`run_covalent_full.py` -> `covalent_full_results.json`)

All six presets × {P1, P1+P2, P4 (vdW-scale sweep λ = 0/0.5/1.0), C3 two-step}:
18/18 runs rc=0.

- P1+P2 NAC discrimination: cys_michael / ser_covalent / lys_targeting /
  tyr_covalent → YES; cys_sn2 (180° backside unsatisfied) and boronic_acid
  (no angle constraint by design) → NO.
- P4 energy monotonicity example (cys_michael): λ=0 → −10.09, λ=0.5 → −5.716,
  λ=1.0 → −4.734 kcal/mol.
- C3 two-step phase-2 pose-set sizes produced for all presets.
- Energy-well scan (cys_michael P1): Gaussian bottom exactly at r₀ = 1.82 Å,
  depth ≈10 kcal/mol; receptor-repulsive component separable from the reactive term.
- `--generate_maps`: **108 AutoGrid 4.2-format affinity maps** written in one
  pass (all AD4 probe types incl. TZ/SQ/MH/JT), zero external dependencies.

## 2. Real Zn-metalloprotein system (4JC, carbonic-anhydrase-like)

Receptor `4JC_rec_zn_fixed.pdbqt` (2556 atoms incl. Zn at −3.38, 0.33, 85.49);
ligand `4JC305_lig.pdbqt` (14 atoms, NA donor). Grid 25 Å³ @ 0.375 Å, exhaust 16,
seed 42, 8 CPUs. NOTE: use the clean receptor (2,083 atoms); a co-crystallized
UNL/UNK-contaminated preparation (1,349 extra atoms) corrupts absolute energies.

| Engine | Score | Best-mode d(Zn–donor) | Metal-coordination report |
|---|---|---|---|
| LKina AD4 + `--metal_mode zn` | −34.49 | **2.13 Å (NA donor)** | `METAL_COORD: Zn donor=NA d=2.13 A`, `METAL_RERANK` present |
| Vina 1.2.7 (vina scoring) | −6.39 | 2.23 Å (OA) | none |
| Vina 1.2.7 (AD4 + LKina-generated maps) | −13.07 | 2.23 Å (OA) | none |

## 3. Reactive covalent presets (6/6, synthetic NAC systems)

Script: `run_reactive_tests.py` -> `reactive_presets_results.json`

| Preset | r0 (Å) | P1 NAC | P1+P2 NAC | P1+P2 d (Å) | P1+P2 angle (°) |
|---|---|---|---|---|---|
| cys_michael | 1.82 | YES | YES | 2.26 | 93.9 |
| cys_sn2 | 1.82 | YES | NO | 2.47 | 108.3 |
| ser_covalent | 1.34 | YES | YES | 2.01 | 99.1 |
| lys_targeting | 1.47 | YES | YES | 2.33 | 94.1 |
| boronic_acid | 1.47 | YES | NO | 2.71 | 110.9 |
| tyr_covalent | 1.38 | YES | YES | 2.03 | 100.8 |

## 4. Backward compatibility (real organic systems, vina scoring)

Re-verified after the engine fixes:

| System | LKina mode-1 score | Vina 1.2.7 mode-1 score | Match |
|---|---|---|---|
| 1HVR / XK2 (HIV-1 protease) | −14.54 | −14.54 | exact |
| 3PTB / BEN (trypsin) | −6.202 | −6.202 | exact |

## 5. Previously reported benchmarks (from LKINA.md design document)

These datasets are part of the engine's published design documentation:
- 292 Zn²⁺ redocking (Santos-Martins protocol): LKina TZ 74.3% top-1 ≤ 2.0 Å
  (AutoDock4Zn 70.5%, standard AD4 58.2%).
- 56 Fe³⁺ (PDBbind v2020, ≤2.5 Å): LKina MH 66.1% vs AD4 48.2%.
- 41 Cu²⁺: LKina JT 63.4% vs AD4 39.0%; JT equatorial ≤0.20 Å 78.3%.
- 28 Cys-Michael covalent: LKina P1+P2 64.3% ≤2.0 Å (Vina 42.9%); NAC 71.4%.

## 6. Supplementary feature tests (run_feature_tests.py)

- **Reactive gradient check (6/6 presets)**: max |analytic − numeric|
  (central difference, eps=1e-4): cys_michael 2.4e-16, cys_sn2 8.0e-09,
  ser_covalent 4.3e-08, lys_targeting 2.9e-08, boronic_acid 9.5e-09,
  tyr_covalent 4.2e-08.
- **Metal Bias (O5)** / **C3 two-step** CLI flags verified.
- Energy-scale sanity (clean 4JC crystal pose, score_only): metal_zn −11.07 =
  metal_soft(w=0.3) −11.07.

## 7. Publication figures (make_figures_v3.py, 2026-08-27)

All figures regenerate from the JSON files in this directory via
`make_figures_v3.py` (repository root of the manuscript folder). Output:
300-dpi PNG **and** vector PDF per figure, Arial, panel labels, journal
sizes (single column 3.5 in, double column 7.2 in).

| Figure | Source JSON | Content |
|---|---|---|
| fig5 | metal_coverage_results_all.json | 110-mode coverage, LKina 110/110 vs Vina 34/110, by-family acceptance |
| fig6 | metal_coverage_results_all.json + feature_family_results.json | coordination-distance histogram (mean 0.20 Å, 108/110 ≤ 0.5 Å); TZ/SQ/MH/JT geometry check |
| fig7 | feature_family_results.json | BVS oxidation-state inference 14/14 |
| fig8 | covalent_full_results.json | four-tier framework 6/6; NAC angle discrimination; Gaussian well scan |
| fig9 | covalent_full_results.json + feature_family_results.json | 108 inline maps; metal-as-ligand QC 0.003 vs 0.903 |
| figS1 | metal_coverage_results_all.json | energy-well landscape over all 110 modes by pseudoatom class (104/110 wells) |
| figS2 | reactive_presets_results.json | P1 vs P1+P2 reactive-distance convergence with r0 reference lines |
