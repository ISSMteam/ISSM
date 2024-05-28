#!/bin/bash
set -eu


# Keeping the following commented line for potential future use.
#git clone https://gitlab.com/adol-c/adol-c.git src

## Environment
#
export CFLAGS="-O2 -L${ISSM_DIR}/externalpackages/petsc/install/lib -lmpi"
export CXXFLAGS="-O2 -std=c++11 -L${ISSM_DIR}/externalpackages/petsc/install/lib -lmpi"

# Cleanup
rm -rf install src
mkdir install src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/ADOL-C.tar.gz" "ADOL-C.tar.gz"

# Unpack source
tar -zxvf ADOL-C.tar.gz

# Configure
cd src
./configure \
	--prefix=${ISSM_DIR}/externalpackages/adolc/install \
	--libdir=${ISSM_DIR}/externalpackages/adolc/install/lib \
	--with-mpi-root=${ISSM_DIR}/externalpackages/petsc/install \
	--enable-ampi \
	--with-ampi=${ISSM_DIR}/externalpackages/adjoinablempi/install \
	--with-soname=adolc \
	--disable-tapedoc-values

# Clean
make clean

# Compile
if [ $# -eq 0 ]; then
	make V=1
else
	make -j $1 V=1
fi

# Install
make V=1 install
