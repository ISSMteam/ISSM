#!/bin/bash
set -eu

#Some cleanup
rm -rf install src build

#Download development version
git clone https://github.com/SciCompKL/adjoint-PETSc.git src
mkdir install
mkdir build

grep -rl 'remove_cvref_t' . | xargs sed -i 's/std::remove_cvref_t<\([A-Za-z0-9_]*\)>/std::remove_cv_t<std::remove_reference_t<\1>>/g'

$ISSM_DIR/externalpackages/cmake/install/bin/cmake src \
    -B build \
    -DBUILD_TESTING=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=$ISSM_DIR/externalpackages/adjointpetsc/install \
    -DCoDiPack_DIR=$ISSM_DIR/externalpackages/codipack/install/cmake \
    -DPETSc_DIR=$ISSM_DIR/externalpackages/petsc/install

if [ $# -eq 0 ]; then
  cmake --build build
else
  cmake --build build -j $1
fi
cmake --install build
