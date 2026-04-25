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
# - Test if -static-libgfortran flag will avoid all of this.
# - Otherwise, refactor this to work with other gfortran installations.
#
echo "Finding libgfortran..."
LIBGFORTRAN=$(find /usr -name libgfortran.* 2>/dev/null | egrep -n libgfortran.a | sed "s/[0-9]*://g" | head -1)
LIBGFORTRAN_ROOT=${LIBGFORTRAN%/*}

GDAL_EXTRA_LIBS="-lm ${HDF5_ROOT}/lib/libhdf5_hl.a ${HDF5_ROOT}/lib/libhdf5.a ${ZLIB_ROOT}/lib/libz.a -ldl -lsqlite3 ${PROJ_ROOT}/lib/libproj.a ${MPI_ROOT}/lib/libmpicxx.a -lstdc++" # See also patch for configuration file ./configs/${VER%.*}/static/cmake/modules/FindGDAL.cmake; for some reason, needed to run `sudo ln -s /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so` after recent upgrades to make it so libstc++ is found
NETCDF_EXTRA_LIBS="${CURL_ROOT}/lib/libcurl.a ${HDF5_ROOT}/lib/libhdf5_hl.a ${HDF5_ROOT}/lib/libhdf5.a ${ZLIB_ROOT}/lib/libz.a -ldl" # See also patch for configuration file ./configs/${VER%.*}/static/cmake/modules/FindNETCDF.cmake)

# Environment
#
export CC=mpicc
export CFLAGS="${CFLAGS} -w"

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
patch ./src/cmake/ConfigUser.cmake < ./configs/${VER%.*}/static/cmake/ConfigUserAdvancedTemplate.cmake.patch
patch ./src/cmake/modules/FindGDAL.cmake < ./configs/${VER%.*}/static/cmake/modules/FindGDAL.cmake.patch
patch ./src/cmake/modules/FindGSHHG.cmake < ./configs/${VER%.*}/static/cmake/modules/FindGSHHG.cmake.patch
patch ./src/cmake/modules/FindNETCDF.cmake < ./configs/${VER%.*}/static/cmake/modules/FindNETCDF.cmake.patch

# Configure
cd src
mkdir build
cd build

# NOTE:
# - The CMake modules used to find and probe the BLAS and LAPACK libraries do
#	not seem to handle the situation where BLAS_LIBRARY and LAPACK_LIBRARY are
#	set but we are working with static libraries
#	(see customized ConfigUser.static.cmake). Using BLAS_LIBRARIES and
#	LAPACK_LIBRARIES is a workaround.
#
cmake \
	-DCMAKE_INSTALL_PREFIX="${PREFIX}" \
	-DBLAS_LIBRARIES="${BLAS_ROOT}/lib/libfblas.a;-L${LIBGFORTRAN_ROOT};-lgfortran" \
	-DCURL_INCLUDE_DIR="${CURL_ROOT}/include" \
	-DCURL_LIBRARY="${CURL_ROOT}/lib/libcurl.a" \
	-DGDAL_EXTRA_LIBS="${GDAL_EXTRA_LIBS}" \
	-DLAPACK_LIBRARIES="${LAPACK_ROOT}/lib/libflapack.a;-L${LIBGFORTRAN_ROOT};-lgfortran" \
	-DNETCDF_EXTRA_LIBS="${NETCDF_EXTRA_LIBS}" \
	..

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
