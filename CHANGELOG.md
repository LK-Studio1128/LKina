# Changelog

All notable changes to LKina relative to AutoDock Vina 1.2.7 are documented here.

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
