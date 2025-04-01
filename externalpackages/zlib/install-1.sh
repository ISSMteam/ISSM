#!/bin/bash
set -eu


# Constants
#
VER="1.2.11"

PREFIX="${ISSM_DIR}/externalpackages/zlib/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/zlib-${VER}.tar.gz" "zlib-${VER}.tar.gz"

# Unpack source
tar -zxvf zlib-${VER}.tar.gz

# Move source into 'src' directory
mv zlib-${VER}/* src/
rm -rf zlib-${VER}

# Configure
cd src
./configure \
	--prefix="${PREFIX}"

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
