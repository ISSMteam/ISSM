#!/bin/bash
set -eu


# NOTE:
# - This configuration depends on the successful installation of METIS via 
# 	$ISSM_DIR/externalpackages/metis/install-5-win-msys2-gcc-msmpi.sh
# - MSMPI_ROOT should be available after etc/environment.sh is 
#	sourced
#

## Constants
#
VER=2.1.0

PREFIX="${ISSM_DIR}/externalpackages/scalapack/install"

# Cleanup
rm -rf ${PREFIX} src
mkdir ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/scalapack-${VER}.tgz" "scalapack-${VER}.tgz"

# Unpack source
tar -zxvf scalapack-${VER}.tgz

# Move source into 'src' directory
mv scalapack-${VER}/* src
rm -rf scalapack-${VER}

# Copy customized source and configuration files to 'src' directory
cp configs/2.1/win/msys2/mingw64/msmpi/CMakeLists.txt src
cp configs/2.1/win/msys2/mingw64/msmpi/BLACS/SRC/src-C.c.in src/BLACS/SRC

# Compile
#
# NOTE:
# - The CMakeLists.txt copied above is really just hacked from the original to 
#	push through library and include paths, especially in the case of the 
#	MS-MPI assets
# - Added -fallow-argument-mismatch option to CMAKE_Fortran_FLAGS in order to 
#	clear "Error: Rank mismatch between actual argument at [...]"
# - Added -fallow-invalid-boz option to CMAKE_Fortran_FLAGS in order to clear 
#	"Error: BOZ literal constant at [...]"
#
# TODO:
# - Modify CMakeLists.txt to allow for more flexibility in configuration
#
cd src
cmake \
	-DCMAKE_INSTALL_PREFIX=${PREFIX} \
	-DBUILD_SHARED_LIBS=ON \
	-DBUILD_STATIC_LIBS=OFF \
	-DCMAKE_C_COMPILER=/mingw64/bin/gcc \
	-DCMAKE_C_FLAGS="-Wno-error=implicit-function-declaration" \
	-DCMAKE_Fortran_COMPILER=/mingw64/bin/gfortran \
	-DCMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz" \
	-DLOCAL_MSMPI=1 \
	-DMSMPI_ROOT="${MSMPI_ROOT}" \
	-DBLAS_LIBRARIES="-L${BLAS_ROOT}/lib -lblas" \
	-DLAPACK_LIBRARIES="-L${LAPACK_ROOT}/lib -llapack" \
	.

make

# Install
mkdir ${PREFIX}/lib
cp lib/libscalapack.* ${PREFIX}/lib

# Create link to shared version of library so that libtool can find it
cd ${PREFIX}/lib
ln -s liblapack.so liblapack.dll
