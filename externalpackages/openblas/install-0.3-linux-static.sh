#!/bin/bash
set -eu


# NOTE: Full LAPACK implementation is included in OpenBLAS
#

## Constants
#
PKG="OpenBLAS"
VER="0.3.7"

## Environment
#
export CC="${ISSM_DIR}/externalpackages/mpich/install/bin/mpicc"
export CXX="${ISSM_DIR}/externalpackages/mpich/install/bin/mpicxx"
export FC="${ISSM_DIR}/externalpackages/mpich/install/bin/mpifort"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/${PKG}-${VER}.tar.gz" "${PKG}-${VER}.tar.gz"

# Unpack source
tar -zxvf ${PKG}-${VER}.tar.gz

# Cleanup
rm -rf build install src
mkdir build install install/lib src

# Move source to 'src' directory
mv ${PKG}-${VER}/* src
rm -rf ${PKG}-${VER}

# Configure
#
cd build
cmake \
	../src

# Compile
make

# Install
cd ..
cp ./build/lib/* ./install/lib
