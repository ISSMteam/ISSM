#!/bin/bash
set -eu


# TODO
# - Move installation of GKlib to $PREFIX to Makefile instead of simply copying 
#	it
#

## Constants
#
VER=5.1.0

PREFIX="${ISSM_DIR}/externalpackages/metis/install"

# Cleanup
rm -rf ${PREFIX} src
mkdir ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/metis-${VER}.tar.gz" "metis-${VER}.tar.gz"

# Unpack source
tar -zxvf metis-$VER.tar.gz

# Move source into 'src' directory
mv metis-${VER}/* src
rm -rf metis-${VER}

# Copy customized source and configuration files to 'src' directory
cp configs/5.1/win/msys2/Makefile src
cp configs/5.1/win/msys2/GKlib/gk_arch.h src/GKlib
cp configs/5.1/win/msys2/GKlib/gk_getopt.h src/GKlib

# Configure
cd src
make config \
	prefix=${PREFIX} \
	shared=1 \
	cc=/mingw64/bin/gcc \
	cxx=/mingw64/bin/g++

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi

# Install GKlib
cp -R GKlib ${PREFIX}

# Create link to lib directory (PETSc, by default, looks for libraries in 
# lib64/ if it detects that 64-bit integers are being used)
#
# NOTE: MSYS2 needs to be run as administrator for this to work.
#
cd ${PREFIX}
ln -s ./lib ./lib64

# Create link to shared version of library so that libtool can find it
cd ${PREFIX}/lib
ln -s libmetis.so libmetis.dll
