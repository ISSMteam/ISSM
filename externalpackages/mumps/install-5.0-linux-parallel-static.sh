#!/bin/bash
set -eu


## NOTE: This install script uses make directly rather than CMake and then make

## Constants
#
PKG="mumps"
VER="5.0.2-p2"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/${PKG}-${VER}.tar.gz" "${PKG}-${VER}.tar.gz"

# Unpack source
tar -zxvf ${PKG}-${VER}.tar.gz

# Cleanup
rm -rf install src
mkdir install install/include install/lib src

# Move source to 'src' directory
mv MUMPS_${VER}/* src
rm -rf MUMPS_${VER}

# Copy customized source and config files to 'src' directory
cp configs/5.0/linux/Makefile.debian.static.PAR src/Makefile.inc

# Compile
cd src
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi

# Install
cd ..
cp src/include/* install/include
cp src/lib/lib*.* install/lib
