#!/bin/bash
set -e


## Constants
#
VER="6.6.0"

PREFIX="${ISSM_DIR}/externalpackages/gmt/install"

# Find certain libraries so we do not have to hardcode them
#
# TODO:
# - Move this to etc/environment.sh
#
echo "Finding libgfortran..."
LIBGFORTRAN=$(find /usr /opt -name libgfortran.* 2>/dev/null | egrep -n libgfortran.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)
LIBGFORTRAN_ROOT=${LIBGFORTRAN%/*}

# Environment
#
export CC=mpicc
export CFLAGS="${CFLAGS} -w"
export LDFLAGS="-Wl,-rpath,'${GDAL_ROOT}/lib' -Wl,-rpath,'${NETCDF_ROOT}/lib'"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://github.com/GenericMappingTools/gmt/archive/refs/tags/${VER}.tar.gz" "gmt-${VER}.tar.gz"

# Unpack source
tar -zxvf gmt-${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Move source to 'src' directory
mv gmt-${VER}/* src
rm -rf gmt-${VER}

# Copy custom configuration files to source
cp ./src/cmake/ConfigUserAdvancedTemplate.cmake ./src/cmake/ConfigUser.cmake

# Patch source
patch ./src/cmake/ConfigUser.cmake < ./configs/${VER%.*}/cmake/ConfigUserAdvancedTemplate.cmake.patch
patch ./src/cmake/modules/ConfigCMake.cmake < ./configs/${VER%.*}/mac/cmake/modules/ConfigCMake.cmake.patch

# Configure
cd src
mkdir build
cd build

# NOTE:
# - There is a CMake variable named CURL_ROOT in src/cmake/ConfigUser.cmake
#	that, ostensibly, allows for supplying the path to curl when it is in a
#	non-standard location. That said, newer versions of CMake will ignore said
#	variable and instead try to find curl itself. Passing in the two options
#	below overrides this behavior.
#
cmake \
	-DCMAKE_INSTALL_PREFIX="${PREFIX}" \
	-DBLAS_LIBRARIES="-lblas;-L${LIBGFORTRAN_ROOT};-lgfortran" \
	-DCURL_INCLUDE_DIR="${CURL_ROOT}/include" \
	-DCURL_LIBRARY="-L${CURL_ROOT}/lib;-lcurl" \
	-DLAPACK_LIBRARIES="-llapack;-L${LIBGFORTRAN_ROOT};-lgfortran" \
	-DLIBGFORTRAN_ROOT="${LIBGFORTRAN_ROOT}" \
	..

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
