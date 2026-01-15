#!/usr/bin/env bash
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export INSTALL_DIR="${PROJECT_ROOT}/install"

micromamba activate delphes-to-hdf5

export PATH="${PROJECT_ROOT}/bin:${INSTALL_DIR}/bin:${PATH}"

export CC="${CONDA_PREFIX}/bin/gcc"
export CXX="${CONDA_PREFIX}/bin/g++"

# add install/lib64/ to LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${INSTALL_DIR}/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
