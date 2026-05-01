#!/usr/bin/env bash
# ============================================================
#  LKina — Linux build script
#
#  LKina is a derivative of AutoDock Vina
#  Original work: Copyright (c) 2006-2010, The Scripps Research Institute
#  Licensed under the Apache License, Version 2.0
#  https://github.com/ccsb-scripps/AutoDock-Vina
#
#  Modifications: metal-aware inline AD4 map generation (LKina extensions)
#  See LICENSE file for full Apache-2.0 license text.
#
#  Requirements: g++ ≥7, libboost-dev ≥1.65
#
#  Debian/Ubuntu:
#    sudo apt-get install g++ libboost-all-dev
#
#  CentOS/RHEL:
#    sudo yum install gcc-c++ boost-devel
#
#  Conda/Miniforge (cross-distro):
#    conda install -c conda-forge boost cxx-compiler
#
#  Usage:
#    chmod +x build_LKina_linux.sh
#    ./build_LKina_linux.sh [--conda-env ENV_NAME]
#
#  Output: build/linux/release/LKina
# ============================================================
set -euo pipefail

CONDA_ENV=""
STATIC=n

# Parse optional --conda-env argument
while [[ $# -gt 0 ]]; do
    case $1 in
        --conda-env) CONDA_ENV="$2"; shift 2 ;;
        --static)    STATIC=y; shift ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

# ---- Locate Boost ------------------------------------------
if [ -n "$CONDA_ENV" ]; then
    # Activate conda environment
    CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
    BASE="$CONDA_PREFIX"
elif [ -n "${CONDA_PREFIX:-}" ]; then
    BASE="$CONDA_PREFIX"
elif [ -d /opt/homebrew ]; then
    BASE=/opt/homebrew
else
    BASE=/usr
fi

BOOST_INCLUDE="$BASE/include"
if [ ! -f "$BOOST_INCLUDE/boost/version.hpp" ]; then
    echo "ERROR: Boost headers not found at $BOOST_INCLUDE" >&2
    echo "Install: sudo apt-get install libboost-all-dev  (or set --conda-env)" >&2
    exit 1
fi

BUILD_DIR="$(dirname "$0")/build/linux/release"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
[ -f dependencies ] || touch dependencies

NCPUS=$(nproc 2>/dev/null || echo 4)
echo "Compiling LKina for Linux ($(uname -m), CPUs=$NCPUS)..."
echo "Boost root: $BASE"

make -j"$NCPUS" LKina \
    BASE="$BASE" \
    BOOST_INCLUDE="$BOOST_INCLUDE" \
    BOOST_STATIC="$STATIC" \
    GPP="g++" \
    C_PLATFORM="-pthread" \
    C_OPTIONS="-O3 -DNDEBUG -std=c++14 -fPIC" \
    BOOST_LIB_VERSION="" \
    2>&1

if [ -f LKina ]; then
    echo ""
    echo "✓  Build successful: $BUILD_DIR/LKina"
    echo "   Version: $(./LKina --version)"
else
    echo "✗  Build failed." >&2
    exit 1
fi
