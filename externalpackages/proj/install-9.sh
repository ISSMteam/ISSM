#!/bin/bash
set -eu


# Constants
#
VER="9.5.0"

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
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://download.osgeo.org/proj/proj-${VER}.tar.gz" "proj-${VER}.tar.gz"

# Unpack source
tar -zxvf proj-${VER}.tar.gz

# Move source into 'src' directory
mv proj-${VER}/* src
rm -rf proj-${VER}

# Configure
cd src
mkdir build
cd build
cmake .. \
	-DCMAKE_INSTALL_PREFIX="${PREFIX}" \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DBUILD_SHARED_LIBS=ON

# Compile and install
cmake --build .
cmake --build . --target install
