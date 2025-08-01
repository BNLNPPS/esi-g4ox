#!/usr/bin/env bash
set -euo pipefail            # Fail fast, fail loud
export DEBIAN_FRONTEND=noninteractive

# ------------------------------------------------------------------
# 1. OS packages
# ------------------------------------------------------------------
apt-get update -y
apt-get install -y \
    libgsl-dev libxpm-dev libxft-dev libtbb-dev binutils cmake dpkg-dev \
    g++ gcc libssl-dev git libx11-dev libxext-dev python3 libgif-dev python3 python3-dev

# ------------------------------------------------------------------
# 2. Build and install ROOT (latest-stable branch)
# ------------------------------------------------------------------
git clone --branch latest-stable --depth 1 https://github.com/root-project/root.git root_src
mkdir -p root_build root_install
cd root_build

cmake -DCMAKE_INSTALL_PREFIX="$(pwd)/../root_install" ../root_src
cmake --build . --target install -j1

# Source ROOT for this shell (optional if you only need it later)
source "$(pwd)/bin/thisroot.sh"
cd ..

# ------------------------------------------------------------------
# 3. Fetch and build RAT-PAC
# ------------------------------------------------------------------
git clone https://github.com/rat-pac/ratpac-setup.git
cd ratpac-setup
./setup.sh -j"$(nproc)"
