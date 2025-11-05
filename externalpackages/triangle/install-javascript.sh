#!/bin/bash
set -eu


# Constants
#
export PREFIX="${ISSM_DIR}/externalpackages/triangle/install" # Set to location where external package should be installed

# Environment
#
export CC=emcc
export CXX=em++
export AR=emar
export RANLIB=emranlib
#export EMCC_DEBUG=1 # Uncomment to enable debugging

# Source Emscripten environment
source ${EMSCRIPTEN_ROOT}/emsdk_env.sh

# Cleanup
rm -rf ${PREFIX}
mkdir -p ${PREFIX} ${PREFIX}/include ${PREFIX}/share src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://github.com/ISSMteam/ExternalPackages/raw/refs/heads/main/triangle.zip" "triangle.zip"

# Unpack source
unzip triangle.zip -d src

# Copy customized source files to 'src' directory
cp configs/makefile src
cp configs/javascript/configure.make src
cp configs/javascript/triangle.h src

# Compile
cd src
make objects

# Install
cd ..
cp src/triangle.o ${PREFIX}/share
cp src/triangle.h ${PREFIX}/include

# Cleanup
rm -rf src
