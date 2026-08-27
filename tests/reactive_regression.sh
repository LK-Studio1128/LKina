#!/usr/bin/env bash
# LKina 17-item automated regression suite (LKINA.md §5.1, groups A–F).
# Self-contained: locates the binary and benchmarks relative to the repo root.
#
# Usage:
#   ./tests/reactive_regression.sh            # auto-locate binary
#   BIN=/path/to/LKina ./tests/reactive_regression.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
BIN="${BIN:-}"
if [[ -z "$BIN" ]]; then
  for c in "$REPO/build/mac/release/LKina" "$REPO/build/linux/release/LKina"; do
    [[ -x "$c" ]] && BIN="$c" && break
  done
fi
if [[ -z "$BIN" ]]; then
  echo "ERROR: LKina binary not found; set BIN=/path/to/LKina" >&2
  exit 2
fi
echo "LKina binary : $BIN"
python3 "$HERE/run_regression.py" --binary "$BIN" --benchmarks "$REPO/benchmarks"
