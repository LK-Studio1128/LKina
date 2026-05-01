#!/usr/bin/env bash
# ============================================================
#  LKina — macOS build script
#
#  LKina is a derivative of AutoDock Vina
#  Original work: Copyright (c) 2006-2010, The Scripps Research Institute
#  Licensed under the Apache License, Version 2.0
#  https://github.com/ccsb-scripps/AutoDock-Vina
#
#  Modifications: metal-aware inline AD4 map generation (LKina extensions)
#  See LICENSE file for full Apache-2.0 license text.
#
#  Requirements: Xcode CLT (clang++), Homebrew, Boost ≥1.71
#
#  Usage:
#    chmod +x build_LKina_mac.sh
#    ./build_LKina_mac.sh
#
#  Output: build/mac/release/LKina
# ============================================================
set -euo pipefail

# ---- Detect Homebrew prefix --------------------------------
if command -v brew >/dev/null 2>&1; then
    BREW_PREFIX=$(brew --prefix)
else
    echo "ERROR: Homebrew not found. Install from https://brew.sh" >&2
    exit 1
fi

# ---- Ensure Boost is installed -----------------------------
if ! brew list boost &>/dev/null; then
    echo "Installing Boost via Homebrew..."
    brew install boost
fi

BUILD_DIR="$(dirname "$0")/build/mac/release"
echo "Build directory: $BUILD_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Copy dependency file if it doesn't exist
[ -f dependencies ] || touch dependencies

echo "Compiling LKina for macOS ($(uname -m))..."
make -j"$(sysctl -n hw.logicalcpu)" LKina \
    BASE="$BREW_PREFIX" \
    BOOST_INCLUDE="$BREW_PREFIX/include" \
    BOOST_STATIC=y \
    GPP="/usr/bin/clang++" \
    C_PLATFORM="-pthread" \
    C_OPTIONS="-O3 -DNDEBUG -std=c++14 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe" \
    2>&1

if [ -f LKina ]; then
    echo ""
    echo "✓  Build successful: $BUILD_DIR/LKina"
    echo "   Version: $(./LKina --version)"
else
    echo "✗  Build failed." >&2
    exit 1
fi
