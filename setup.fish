#!/usr/bin/env fish
set -xg PROJECT_ROOT (dirname (readlink -m (status --current-filename)))
set -xg INSTALL_DIR {$PROJECT_ROOT}/install

micromamba activate delphes-to-hdf5

fish_add_path -g {$PROJECT_ROOT}/bin
fish_add_path -g {$INSTALL_DIR}/bin

set -xg CC {$CONDA_PREFIX}/bin/gcc
set -xg CXX {$CONDA_PREFIX}/bin/g++

# add install/lib64/ to LD_LIBRARY_PATH
set -xga LD_LIBRARY_PATH {$INSTALL_DIR}/lib64
