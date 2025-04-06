#!/bin/bash
set -eu


## Constants
#
VER="3.31.6"

PREFIX="${ISSM_DIR}/externalpackages/cmake/install" # Set to location where external package should be installed

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://cmake.org/files/v${VER%.*}/cmake-${VER}.tar.gz" "cmake-${VER}.tar.gz"

# Unpack source
tar -zxvf cmake-${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX}

# Move source into 'src' directory
mv cmake-${VER} src

# Configure
cd src
#./bootstrap \
#	--prefix=${PREFIX} # Breaks on ronne
./configure \
	--prefix="${PREFIX}"

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make install -j $1
fi
