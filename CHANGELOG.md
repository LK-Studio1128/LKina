# Changelog

All notable changes to LKina relative to AutoDock Vina 1.2.7 are documented here.

## LKina 1.0.2

### Fixed
- **Metal-mode TZ well 6× too deep (major scoring bug)**: `nbp_r_eps`
  override epsilons were applied unweighted. Real AutoGrid multiplies
  `epsij *= AD4.coeff_vdW` when building maps (mainpost1.28.cpp:1786), so the
  official `nbp_r_eps 0.25 23.2135 12 6 NA TZ` override yields an actual well
  depth of **−3.86 kcal/mol**, not −23.21. LKina now applies the same
  weighting, consistent with the standard-vdW branch.
- **TZ/SQ/MH/JT pseudoatom parameters**: built-in parameters wrongly carried
  the unweighted override eps (`Rii=0.25, epsii=23.2/20/18/10`), giving every
  non-override probe a spurious LJ well against injected pseudoatoms. The
  official AD4Zn.dat TZ entry is `Rii=1.0, epsii=0.0` (all contacts flow
  through explicit nbp overrides); LKina now uses the official neutral
  parameters.
- **`AG4_EINTCLAMP` semantics restored to AutoGrid**: the official generator
  only clamps *pairwise LJ table values from above* at 100000
  (`autocomm.h:89`) and never clamps the summed affinity-map value (official
  zinc-example maps reach +200437). LKina's ±1000 two-sided clamp (v1.0.1) is
  reverted accordingly. Note: this supersedes the v1.0.1 clamp change, which
  was based on a misreading of AutoGrid's `MAXVALUE` behaviour.

### Verified
- 3HS4 (CA II + AZM) crystal pose, same box: std AD4 −7.453 (unchanged),
  zn mode −23.914 → **−11.324 kcal/mol** (experimental affinity ≈−10.8).
- 3HS4 global docking best: −35.742 → **−13.319**; Zn–ligand 2.08 Å.
- 4JC global docking best: −34.858 → **−11.478**; Table-6 rerun −11.514.
- Non-metal backward compatibility: 1HVR −14.54 / 3PTB −6.202 kcal/mol
  reproduce Vina 1.2.7 exactly; 6OIM covalent battery numerically unchanged.
- Full parameter-sweep battery (exhaustiveness/seed/soft-weight/metal-bias/
  reactive-strength) re-run on the fixed binary; seed spread ≤0.010 kcal/mol.

## LKina 1.0.1

### Fixed
- **AG4 affinity-map clamp (energy-spike fix)**: `AG4_EINTCLAMP` lowered from
  ±100000 to **±1000 kcal/mol** and a final per-grid-point clamp added in
  `ag4_compute_maps`, exactly matching AutoGrid 4.2's `MAXVALUE` behaviour.
  Previously, unphysical repulsive spikes (e.g. a TZ pseudoatom sampled near
  r = 0 on a 0.375 Å grid, measured +14630 kcal/mol) could propagate into the
  affinity maps and dominate poses.
- **Signed-clamp in the LJ pair table**: `ag4_build_lj_table` now clamps the
  lower side of the pair energy as well (`max(-CLAMP, min(CLAMP, E))`),
  preventing deep-well undershoot in addition to the previous upper clamp.

### Added
- **`--no_auto_metal` flag**: opts out of automatic metal-mode detection from
  receptor/ligand PDBQT atom types, keeping pure AD4 scoring. Recommended when
  metal-mode absolute energies (different scale than standard AD4/Vina) are not
  wanted. A stderr notice is printed whenever auto-detection fires.

### Verified
- Non-metal backward compatibility: BEN benchmark score-only reproduces
  AutoDock Vina numerics (−5.806 kcal/mol, identical decomposition).
- Metal path smoke test: 4JC Zn-shell receptor, `--metal_mode zn
  --generate_maps` writes 108 probes × 64³ grids and completes scoring.
- Full 110-mode synthetic benchmark re-run on the fixed binary
  (LKina 110/110, mean |d−d₀| = 0.20 Å; see manuscript benchmarks).

## LKina 1.0.0 (initial release)

### Added
- `--generate_maps` flag: inline AD4 affinity-map generation without an
  external `autogrid4` binary
- `ag4_engine.cpp` / `ag4_engine.h`: independent C++ reimplementation of
  the AutoGrid 4.2 pairwise interaction scoring function
- `embedded_ad4_grid.cpp` / `embedded_ad4_grid.h`: drives `ag4_engine` to
  write AutoGrid4-compatible `.map` files into a temporary directory, then
  passes them to the existing Vina AD4 cache loader
- `ad4_parameter_data.cpp` / `ad4_parameter_data.h`: embedded AD4 force
  field parameters extended with 70+ metal atom types (transition metals,
  lanthanides, metalloids: As, B, Ge, Sb, Se, Te, Tg, …)
- Build scripts: `build_LKina_mac.sh`, `build_LKina_linux.sh`,
  `build_LKina_win.bat` with Apache-2.0 attribution headers
- `NOTICE` file documenting upstream AutoDock Vina attribution

### Changed
- `src/main/main.cpp`: added `--generate_maps` option; extended metal
  atom-type whitelist
- `build/makefile_common`: added `LKina` target alongside `vina` and
  `vina_split`
- `ad4cache.cpp`: extended to accept maps produced by the inline engine

### Unchanged from AutoDock Vina 1.2.7
- Vina scoring function
- Vinardo scoring function
- Standard AD4 scoring (external map files)
- Multithreading, gradient optimization, conformational search
- All CLI flags except the new `--generate_maps`
