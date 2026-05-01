# LKina

[![Release](https://img.shields.io/github/v/release/LK-Studio1128/LKina)](https://github.com/LK-Studio1128/LKina/releases/latest)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CI](https://github.com/LK-Studio1128/LKina/actions/workflows/compile-binaries.yaml/badge.svg)](https://github.com/LK-Studio1128/LKina/actions/workflows/compile-binaries.yaml)

**Repository:** <https://github.com/LK-Studio1128/LKina>

**LKina** is a metal-aware molecular docking engine derived from
[AutoDock Vina 1.2.7](https://github.com/ccsb-scripps/AutoDock-Vina).  
It extends Vina with **inline AD4 affinity-map generation**, enabling
first-class docking of metal-containing ligands (Zn, Fe, Mg, Ca, Cu, Mn,
Au, Pt, and 60+ additional metal atom types) without requiring an external
`autogrid4` binary.

## Key Features

- **Inline AD4 map generation** — pass `--scoring ad4 --generate_maps` and
  LKina computes receptor affinity grids on-the-fly using a built-in
  reimplementation of the AutoGrid 4.2 scoring function
- **70+ metal atom types** — full `_LKINA_SUPPORTED_METALS` whitelist
  covering transition metals, lanthanides, metalloids (As, B, Ge, Sb, Se,
  Te, Tg, …)
- **Drop-in Vina replacement** — same CLI flags; all existing Vina scoring
  functions (Vina, Vinardo, AD4) remain available
- **No external autogrid4 dependency** — maps are generated internally at
  docking time
- Builds on **macOS (arm64/x86_64)**, **Linux (x86_64)**, and
  **Windows (MSYS2/MinGW-w64)**

## Build

### macOS
```bash
chmod +x build_LKina_mac.sh
./build_LKina_mac.sh
# Output: build/mac/release/LKina
```

### Linux
```bash
chmod +x build_LKina_linux.sh
./build_LKina_linux.sh
# Output: build/linux/release/LKina
```

### Windows (MSYS2 required)
```
build_LKina_win.bat
# Output: build/win/release/LKina.exe
```

**Dependencies:** g++/clang++ ≥ 7, Boost ≥ 1.65 (program_options,
thread, serialization, filesystem)

## Usage

```bash
# Standard Vina scoring (unchanged)
LKina --receptor rec.pdbqt --ligand lig.pdbqt \
      --center_x 0 --center_y 0 --center_z 0 \
      --size_x 20 --size_y 20 --size_z 20 --out out.pdbqt

# Inline AD4 metal docking (LKina extension)
LKina --scoring ad4 --generate_maps \
      --receptor rec.pdbqt --ligand metal_lig.pdbqt \
      --center_x 14.68 --center_y 32.38 --center_z 10.64 \
      --size_x 25 --size_y 25 --size_z 25 --out out.pdbqt
```

## Documentation

- [LKINA.md](LKINA.md) — full design paper (English abstract + Chinese body) covering the
  113-type atom system, TZ/SQ/MH/JT pseudoatoms, BVS oxidation-state inference,
  reactive covalent framework P1–P4, C3 two-step strategy, 6 reaction presets,
  Metal Bias (O5) and benchmark results. Quick-start and cheat sheets are in
  Appendices A–F.
- [CHANGELOG.md](CHANGELOG.md) — release notes vs AutoDock Vina 1.2.7.
- [RELEASE.md](RELEASE.md) — packaging and GitHub Release process.

## Contributing

Bug reports and pull requests are welcome at <https://github.com/LK-Studio1128/LKina>.
See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

LKina uses a **dual-license** structure:

| Component | License | Files |
|-----------|---------|-------|
| AutoDock Vina core | Apache-2.0 | all files **except** those listed below |
| AG4 engine & parameter data | **GPL-3.0-or-later** | `ag4_engine.*`, `embedded_ad4_grid.*`, `ad4_parameter_data.*` |

Because these components are linked into a single binary, the **LKina
executable as a whole is distributed under GPL-3.0-or-later**.

- [LICENSE](LICENSE) — GPL-3.0 full text (primary distribution license)
- [LICENSE.Apache-2.0](LICENSE.Apache-2.0) — Apache-2.0 full text (Vina-origin files)
- [COPYING](COPYING) — dual-license scope, file map and third-party inventory
- [NOTICE](NOTICE) — upstream copyright, attribution and modification log

The AG4 files reference AutoGrid4 (GPL-2.0 or later,
© The Scripps Research Institute). Apache-2.0 is compatible with
GPL-3.0, so the combined binary can legally be distributed under GPL-3.0.

## Citations

If you use LKina in published work, please cite the upstream Vina papers:

- J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli. (2021).
  AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and
  Python Bindings. *J. Chem. Inf. Model.*
  <https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203>
- O. Trott and A. J. Olson. (2010). AutoDock Vina: improving the speed and
  accuracy of docking with a new scoring function, efficient optimization,
  and multithreading. *J. Comput. Chem.* 31(2), 455–461.
  <https://onlinelibrary.wiley.com/doi/10.1002/jcc.21334>
- G. Morris et al. (2009). AutoDock4 and AutoDockTools4: Automated docking
  with selective receptor flexibility. *J. Comput. Chem.* 30, 2785–2791.
