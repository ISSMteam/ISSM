#!/bin/bash
set -eu


## Constants
#
VER="2.7"

PREFIX="${ISSM_DIR}/externalpackages/gsl/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/gsl-${VER}.tar.gz" "gsl-${VER}.tar.gz"

# Unpack source
tar -zxvf gsl-${VER}.tar.gz

# Move source to 'src' directory
mv gsl-${VER}/* src
rm -rf gsl-${VER}

# Configure
cd src
./configure \
	--prefix="${PREFIX}" \
	--disable-shared \
	--with-pic

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
