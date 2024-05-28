#!/bin/bash
set -eu


# TODO:
# - Modify configuration/makefiles so that dll rather than so is produced.
#

## Constants
#
VER="3.9.0"

PREFIX="${ISSM_DIR}/externalpackages/lapack/install"

# Cleanup
rm -rf ${PREFIX} build src
mkdir ${PREFIX} build src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/lapack-${VER}.tar.gz" "lapack-${VER}.tar.gz"

# Unpack source
tar -zxvf lapack-${VER}.tar.gz

# Move source into 'src' directory
mv lapack-${VER}/* src
rm -rf lapack-${VER}

# Configure
#
cd build
cmake \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_C_COMPILER=/mingw64/bin/gcc \
	-DCMAKE_Fortran_COMPILER=/mingw64/bin/gfortran \
	-DBLAS_LIBRARIES="-L${BLAS_ROOT}/lib -lblas" \
	../src

# Compile
make

# Copy libraries to lib directory (on Windows, CMake installs .dll files next
# to the .exe files that need them, making tests easy to run, so let's not 
# change the CMake configuration).
mkdir ${PREFIX}/lib
cp lib/liblapack.* ${PREFIX}/lib

# # Create link to versioned library
# cd ${PREFIX}/lib
# ln -s msys-lapack-3.dll msys-lapack.dll

# Create link to shared version of library so that libtool can find it
cd ${PREFIX}/lib
ln -s liblapack.so liblapack.dll
