# Contributing to LKina

Thank you for your interest in contributing to LKina!

## Bug Reports

Please open a GitHub Issue at <https://github.com/LK-Studio1128/LKina/issues> with:

- LKina version (`LKina --version`)
- Operating system and architecture
- Minimal command line that reproduces the issue
- Full error output

## Pull Requests

1. Fork the repository and create a feature branch from `main`.
2. Ensure your changes build cleanly on at least one supported platform:
   ```bash
   # macOS
   ./build_LKina_mac.sh
   # Linux
   ./build_LKina_linux.sh
   ```
3. If you add or modify the AG4 engine or parameter files, make sure the license
   header in the affected file(s) says **GPL-3.0-or-later** (not Apache-2.0).
4. Submit a pull request with a clear description of the change and the
   motivation.

## Code Style

LKina follows the existing AutoDock Vina C++ code style:

- C++14 standard (`-std=c++14`)
- 4-space indentation
- Descriptive variable names; no single-letter temporaries in non-trivial loops

## License

By submitting a pull request you agree that your contribution will be licensed
under the same terms as the file(s) you modify:

- Modifications to the AG4 engine / parameter files → **GPL-3.0-or-later**
- Modifications to other Vina-derived files → **Apache-2.0**
- New files must include an appropriate SPDX license header

See [COPYING](COPYING) and [LICENSE](LICENSE) for details.
