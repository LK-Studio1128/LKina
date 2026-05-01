# Changelog

All notable changes to LKina relative to AutoDock Vina 1.2.7 are documented here.

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
