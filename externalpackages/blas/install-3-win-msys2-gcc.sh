#!/bin/bash
set -eu


## Constants
#
VER=3.8.0

PREFIX="${ISSM_DIR}/externalpackages/blas/install"

MODULE="blas"
IMP_LIB_NAME="lib${MODULE}.dll.a"
LIB_NAME="msys-${MODULE}.dll"

# Cleanup
rm -rf ${PREFIX} src
mkdir ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/blas-${VER}.tgz" "blas-${VER}.tgz"

# Unpack source
tar -zxvf blas-${VER}.tgz

# Move source into 'src' directory
mv BLAS-${VER}/* src
rm -rf BLAS-${VER}

# Copy customized source and configuration files to 'src' directory
cp configs/3.8/win/msys2/gcc/make.inc src
cp configs/3.8/win/msys2/gcc/Makefile src

# Compile
cd src
make

# Install
mkdir ${PREFIX}/lib
cp ${IMP_LIB_NAME} ${PREFIX}/lib
cp ${LIB_NAME} ${PREFIX}/lib

# Create link to shared version of library so that libtool can find it
cd ${PREFIX}/lib
ln -s ./${LIB_NAME} ./libblas.dll
