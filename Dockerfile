# LKina — reproducible build & run container
#
# Build:   docker build -t lkina .
# Run:     docker run --rm lkina --help
#          docker run --rm -v "$PWD:/data" lkina --receptor /data/rec.pdbqt \
#                     --ligand /data/lig.pdbqt --metal_mode zn --out /data/out.pdbqt
# Test:    docker run --rm --entrypoint /src/tests/reactive_regression.sh lkina
#
# LKina is a derivative of AutoDock Vina (Apache-2.0, The Scripps Research
# Institute); LKina extensions are GPL-3.0-or-later. This container installs
# only system libraries (boost headers, libgomp) and builds from the local
# source tree.

# ---------- stage 1: build ----------
FROM debian:bookworm-slim AS builder

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        g++ make ca-certificates libboost-all-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

# g++ >= 7 and Boost >= 1.65 headers are the only requirements
# (see build_LKina_linux.sh). -fopenmp is added by the project makefile.
RUN ./build_LKina_linux.sh

# smoke check inside the builder
RUN ./build/linux/release/LKina --version

# ---------- stage 2: runtime ----------
FROM debian:bookworm-slim

LABEL org.opencontainers.image.title="LKina" \
      org.opencontainers.image.description="Metal-aware and covalent-reactive molecular docking engine extending AutoDock Vina 1.2.7" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.source="https://github.com/LK-Studio1128/LKina"

# OpenMP runtime is the only shared-library dependency of the binary
RUN apt-get update \
 && apt-get install -y --no-install-recommends libgomp1 ca-certificates \
 && rm -rf /var/lib/apt/lists/* \
 && useradd -m -u 1000 lkina

COPY --from=builder /src/build/linux/release/LKina /usr/local/bin/LKina

USER lkina
WORKDIR /data

ENTRYPOINT ["LKina"]
CMD ["--help"]
