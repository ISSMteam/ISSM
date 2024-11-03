#!/bin/bash
set -eu


## Constants
#
VER="2.7"

PREFIX="${ISSM_DIR}/externalpackages/gsl/install" # Set to location where external package should be installed

## Environment
#
export CC=emcc
export CXX=em++
export AR=emar
export RANLIB=emranlib
#export EMCC_DEBUG=1 # Uncomment to enable debugging

# Source Emscripten environment
source ${EMSCRIPTEN_ROOT}/emsdk_env.sh

# Issue with variadic function signatures.
#export CFLAGS=-DSTDC_HEADERS

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://ftp.gnu.org/gnu/gsl/gsl-${VER}.tar.gz" "gsl-${VER}.tar.gz"

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

#Compile gsl
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
