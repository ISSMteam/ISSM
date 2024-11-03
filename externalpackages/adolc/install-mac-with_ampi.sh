#!/bin/bash
set -eu


## Constants
#
VER="2.7.2"

## Environment
#
export CFLAGS="-O2 -L${ISSM_DIR}/externalpackages/petsc/install/lib -lmpi"
export CXXFLAGS="-O2 -L${ISSM_DIR}/externalpackages/petsc/install/lib -lmpi"

# Cleanup
rm -rf install src
mkdir install src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://github.com/coin-or/ADOL-C/archive/refs/tags/releases/${VER}.tar.gz" "ADOL-C.tar.gz"

# Unpack source
tar -zxvf ADOL-C.tar.gz

# Move source to 'src' directory
mv ADOL-C-releases-${VER}/* src/
rm -rf ADOL-C-releases-${VER}

# Configure
cd src
./configure \
	--prefix="${ISSM_DIR}/externalpackages/adolc/install" \
	--libdir="${ISSM_DIR}/externalpackages/adolc/install/lib" \
	--with-mpi-root="${ISSM_DIR}/externalpackages/petsc/install" \
	--enable-ampi \
	--with-ampi="${ISSM_DIR}/externalpackages/adjoinablempi/install" \
	--with-soname=adolc \
	--disable-tapedoc-values

# Compile and install
if [ $# -eq 0 ]; then
	make V=1
	make V=1 install
else
	make V=1 -j $1
	make V=1 install -j $1
fi
