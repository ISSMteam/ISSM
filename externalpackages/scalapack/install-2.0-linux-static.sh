#!/bin/bash
set -eu


## NOTE: This install script uses make directly rather than CMake and then make

## Constants
#
PKG="scalapack"
VER="2.0.2"

## Environment
#
export MPI_BASE_DIR="${ISSM_DIR}/externalpackages/mpich/install"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/${PKG}-${VER}.tgz" "${PKG}-${VER}.tgz"

# Unpack source
tar -zxvf ${PKG}-${VER}.tgz

# Cleanup
rm -rf build install src
mkdir build install install/lib src

# Move source to 'src' directory
mv ${PKG}-${VER}/* src
rm -rf ${PKG}-${VER}

# Copy customized source and config files to 'src' directory
cp configs/2.0/linux/SLmake.inc.static src/SLmake.inc

# Compile
cd src
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi

# Install
cd ..
cp ./src/lib*.* ./install/lib
