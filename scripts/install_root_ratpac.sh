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

apt-get install -y libgl1-mesa-dev libglu1-mesa-dev
apt-get install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools

# ------------------------------------------------------------------
# 2. Build and install ROOT (latest-stable branch)
# ------------------------------------------------------------------

spack install root

# ------------------------------------------------------------------
# 3. Fetch and build RAT-PAC
# ------------------------------------------------------------------
git clone https://github.com/rat-pac/ratpac-setup.git
cd ratpac-setup
./setup.sh -j"$(nproc)"

./ratpac/build/bin/rat -o output.root -l log.txt ./ratpac/install/share/RAT/macros/examples/electron.mac
