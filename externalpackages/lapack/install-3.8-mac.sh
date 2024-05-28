#!/bin/bash
set -eu


# NOTE: This installation script will build both BLAS and LAPACK libraries
#

## Constants
#
VER="3.8.0"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/lapack-${VER}.tar.gz" "lapack-${VER}.tar.gz"

# Unpack source
tar -zxvf lapack-$VER.tar.gz

# Cleanup
rm -rf build install src
mkdir build install install/lib src

# Move source to 'src' directory
mv lapack-$VER/* src
rm -rf lapack-$VER

# Copy customized configuration files to 'src' directory
cp configs/mac/3.8/CMakeLists.txt src/CMakeLists.txt

# Configure
#
cd build
cmake \
	-DBUILD_SHARED_LIBS=ON \
	../src

# Compile
make

# Install
cd ..
cp ./build/lib/* ./install/lib
