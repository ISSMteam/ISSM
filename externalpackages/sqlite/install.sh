#!/bin/bash
set -eu


# Constants
#
VER="3.51.2"

PREFIX="${ISSM_DIR}/externalpackages/sqlite/install" # Set to location where external package should be installed

# Environment
#
export CFLAGS="-DSQLITE_ENABLE_COLUMN_METADATA=1"

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://github.com/sqlite/sqlite/archive/refs/tags/version-${VER}.tar.gz" "sqlite-version-${VER}.tar.gz"

# Unpack source
tar -zxvf sqlite-version-${VER}.tar.gz

# Move source into 'src' directory
mv sqlite-version-${VER}/* src
rm -rf sqlite-version-${VER}

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
