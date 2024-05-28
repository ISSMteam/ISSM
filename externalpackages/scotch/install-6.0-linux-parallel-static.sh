#!/bin/bash
set -eu


## NOTE: This install script uses make directly rather than CMake and then make

## Constants
#
PKG="scotch"
VER="6.0.9"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/${PKG}_${VER}.tar.gz" "${PKG}_${VER}.tar.gz"

# Unpack source
tar -zxvf ${PKG}_${VER}.tar.gz

# Cleanup
rm -rf install src
mkdir install install/include install/lib src

# Move source to 'src' directory
mv ${PKG}_${VER}/* src
rm -rf ${PKG}_${VER}

# Copy customized source and config files to 'src' directory
cp configs/6.0/linux/Makefile.inc.linux-parallel-static src/src/Makefile.inc

# Compile
cd src/src
if [ $# -eq 0 ]; then
	make ptscotch
	make ptesmumps
else
	make -j $1 ptscotch
	make -j $1 ptesmumps
fi

# Install
cd ../..
cp src/include/* install/include
cp src/lib/lib*.* install/lib
