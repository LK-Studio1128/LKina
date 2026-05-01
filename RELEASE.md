# LKina Release Process

## Pre-release checklist

- [ ] Update `LKINA_VERSION` in `build/makefile_common`
- [ ] Update `CHANGELOG.md` with all notable changes
- [ ] Verify build on all three platforms (macOS arm64, Linux x86-64, Windows x64)
- [ ] Run `test_LKina_metals.sh` and confirm it passes
- [ ] Confirm `--version` output matches the new version string

## Building release binaries

### macOS arm64

```bash
./build_LKina_mac.sh
# → build/mac/release/LKina
```

### Linux x86-64

```bash
./build_LKina_linux.sh
# → build/linux/release/LKina
```

### Windows x64 (MSYS2 MinGW-w64)

```
build_LKina_win.bat
# → build/win/release/LKina.exe
```

## Packaging

Each release asset must bundle the binary together with:

- `LICENSE` (GPL-3.0 full text — primary distribution license)
- `LICENSE.Apache-2.0` (Apache-2.0 full text — Vina-origin files)
- `COPYING` (dual-license scope and file map)
- `NOTICE` (upstream attribution and modification log)
- `README.md`
- `CHANGELOG.md`

```bash
# Example: macOS arm64 package
mkdir LKina-macos-arm64
cp build/mac/release/LKina \
   LICENSE LICENSE.Apache-2.0 COPYING NOTICE \
   README.md CHANGELOG.md \
   LKina-macos-arm64/
tar -czf LKina-macos-arm64.tar.gz LKina-macos-arm64/
shasum -a 256 LKina-macos-arm64.tar.gz
```

Repeat for Linux (`.tar.gz`) and Windows (`.zip`).

## GitHub Release

1. Create and push the version tag:

   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```

2. Create a GitHub Release at <https://github.com/LK-Studio1128/LKina/releases/new>:
   - **Tag**: `v1.0.0`
   - **Title**: `LKina v1.0.0`
   - **Body**: paste relevant section from `CHANGELOG.md`
   - **Assets**: upload the three platform archives and `SHA256SUMS.txt`

## License compliance notes

LKina is a **derivative work** of AutoDock Vina (Apache-2.0) and AutoGrid4
(GPL-2.0-or-later). Because GPL-3.0-licensed AG4 engine files are linked into
the same binary as Apache-2.0 Vina files, the combined LKina executable must
be distributed under **GPL-3.0**.

- Apache-2.0 is compatible with GPL-3.0 (FSF guidance).
- Every release archive must include `LICENSE`, `COPYING`, and `NOTICE`.
- Source code must remain publicly available at the GitHub repository.
