#!/bin/bash
set -eu


## Constants
#
VER="2.7.2"

# Cleanup
rm -rf install src

# Keeping the following commented line for potential future use.
#git clone https://gitlab.com/adol-c/adol-c.git src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://github.com/coin-or/ADOL-C/archive/refs/tags/releases/${VER}.tar.gz" "ADOL-C.tar.gz"
# Unpack source
tar -zxvf ADOL-C.tar.gz

# Move source to 'src' directory
mv ADOL-C-releases-${VER}/* src/
rm -rf ADOL-C-releases-${VER}

# Configure and compile
cd src
./configure \
	--prefix=$ISSM_DIR/externalpackages/adolc/install  \
	--libdir=$ISSM_DIR/externalpackages/adolc/install/lib \
	--with-mpi-root=$ISSM_DIR/externalpackages/mpich/install \
	--enable-ampi \
	--with-ampi=$ISSM_DIR/externalpackages/adjoinablempi/install \
	--with-soname=adolc \
	--disable-tapedoc-values

make clean
if [ $# -eq 0 ]; then
	make V=1
else
	make -j $1 V=1
fi
make V=1 install
