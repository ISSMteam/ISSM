#!/bin/bash
set -eu


# NOTE:
# - This configuration depends on the successful installation of METIS via 
# 	$ISSM_DIR/externalpackages/metis/install-5-win-msys2-gcc-msmpi.sh
# - METIS_ROOT and MSMPI_ROOT should be available after etc/environment.sh is 
#	sourced
#

## Constants
#
VER=4.0.3

PREFIX="${ISSM_DIR}/externalpackages/parmetis/install"

GKLIB_ROOT="${METIS_ROOT}/GKlib"

# Cleanup
rm -rf ${PREFIX} src
mkdir ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/parmetis-${VER}.tar.gz" "parmetis-${VER}.tar.gz"

# Unpack source
tar -zxvf parmetis-${VER}.tar.gz

# Move source into 'src' directory
mv parmetis-${VER}/* src
rm -rf parmetis-${VER}

# Copy customized source and configuration files to 'src' directory
cp configs/4.0/win/msys2/CMakeLists.txt src
cp configs/4.0/win/msys2/Makefile src
cp configs/4.0/win/msys2/include/parmetis.h src/include
cp configs/4.0/win/msys2/mingw64/libparmetis/CMakeLists.txt src/libparmetis

# Configure
cd src
make config \
	prefix=${PREFIX} \
	shared=1 \
	cc=/mingw64/bin/gcc \
	cxx=/mingw64/bin/g++ \
	metis_path=${METIS_ROOT} \
	gklib_path=${GKLIB_ROOT} \
	local_msmpi=1 \
	msmpi_root=${MSMPI_ROOT}

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi

# Create link to lib directory (PETSc, by default, looks for libraries in 
# lib64/ if it detects that 64-bit integers are being used)
cd ${PREFIX}
ln -s ./lib ./lib64

# Create link to shared version of library so that libtool can find it
cd ${PREFIX}/lib
ln -s libparmetis.so libparmetis.dll
