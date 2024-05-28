#!/bin/bash
set -eu


# Constants
#
VER="3300100"

PREFIX="${ISSM_DIR}/externalpackages/sqlite/install" # Set to location where external package should be installed

# Environment
#
export CFLAGS="-DSQLITE_ENABLE_COLUMN_METADATA=1"

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/sqlite-autoconf-${VER}.tar.gz" "sqlite-autoconf-${VER}.tar.gz"

# Unpack source
tar -zxvf sqlite-autoconf-${VER}.tar.gz

# Move source into 'src' directory
mv sqlite-autoconf-${VER}/* src
rm -rf sqlite-autoconf-${VER}

# Configure
cd src
./configure \
	--prefix="${PREFIX}" \
	--enable-fast-install \
	--enable-static=no

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
