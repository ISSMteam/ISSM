#!/bin/bash
set -eu


# Constants
#
VER="8.1.0"

PREFIX="${ISSM_DIR}/externalpackages/proj/install" # Set to location where external package should be installed

## Environment
#
export CC=mpicc
export CXX=mpicxx

# NOTE: On macOS, SQLite3 should be installed by default, but PROJ currently 
# requires,
#
#	SQLITE3_LIBS="-lsqlite3".
#
# On Ubuntu Linux, install the SQLite3 binary, libraries and headers with,
#
#	apt-get install sqlite3 libsqlite3-dev
#
export SQLITE3_LIBS="-lsqlite3"

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/proj-${VER}.tar.gz" "proj-${VER}.tar.gz"

# Unpack source
tar -zxvf proj-${VER}.tar.gz

# Move source into 'src' directory
mv proj-${VER}/* src
rm -rf proj-${VER}

# Configure
cd src
./configure \
	--prefix="${PREFIX}" \
	--disable-dependency-tracking \
	--enable-fast-install \
	--disable-shared \
	--disable-tiff

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
