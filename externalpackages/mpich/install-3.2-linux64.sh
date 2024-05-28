#!/bin/bash
set -eu


## Constants
#
VER="3.2"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/mpich-${VER}.tar.gz" "mpich-${VER}.tar.gz"

# Unpack source
tar -zxvf  mpich-$VER.tar.gz

# Cleanup
rm -rf src install
mkdir src install

# Move source into 'src' directory
mv mpich-$VER/* src
rm -rf mpich-$VER

# Configure
cd src
./configure \
	--prefix="${ISSM_DIR}/externalpackages/mpich/install" \
	--enable-shared

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi

# Return to initial directory
cd ..
