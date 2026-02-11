#!/bin/bash
set -e

# NOTE: After GMT version 6.0.0, we must build a shared copy of libgmt for our 
# distributables as GMT moved from individual executables to a single 
# executable with modules. For example,
#
#	gmtselect -> gmt select
#
# These modules are not compiled into the static version of libgmt (and likely) 
# never will be.
#

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
LIBGFORTRAN=$(find /usr /opt -name libgfortran.* 2>/dev/null | egrep -n libgfortran.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)
LIBGFORTRAN_ROOT=${LIBGFORTRAN%/*}
echo "Finding libgcc..."
LIBGCC=$(find ${LIBGFORTRAN_ROOT} -name libgcc.* 2>/dev/null | egrep -n libgcc.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)
echo "Finding libglib..."
LIBGLIB=$(find /usr /opt -name libglib-2.0.* 2>/dev/null | egrep -n libglib-2.0.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)
LIBGLIB_ROOT=${LIBGLIB%/*}
GLIB_H=$(find /usr /opt -name glib.h 2>/dev/null | head -1)
GLIB_H_DIR=${GLIB_H%/*}
GLIBCONFIG_H=$(find /opt -name glibconfig.h 2>/dev/null | head -1)
GLIBCONFIG_H_DIR=${GLIBCONFIG_H%/*}
echo "Finding libintl..."
LIBINTL=$(find /usr /opt -name libintl.* 2>/dev/null | egrep -n libintl.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)
echo "Finding libpcre..."
LIBPCRE=$(find /usr /opt -name libpcre2-8* 2>/dev/null | egrep -n libpcre2-8.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)

GDAL_EXTRA_LIBS="-lm ${HDF5_ROOT}/lib/libhdf5_hl.a ${HDF5_ROOT}/lib/libhdf5.a ${ZLIB_ROOT}/lib/libz.a -ldl -lsqlite3 ${PROJ_ROOT}/lib/libproj.a ${MPI_ROOT}/lib/libmpicxx.a -lc++" # See also patch for configuration file ./configs/${VER%.*}/static/cmake/modules/FindGDAL.cmake
NETCDF_EXTRA_LIBS="${CURL_ROOT}/lib/libcurl.a ${HDF5_ROOT}/lib/libhdf5_hl.a ${HDF5_ROOT}/lib/libhdf5.a ${ZLIB_ROOT}/lib/libz.a" # See also patch for configuration file ./configs/${VER%.*}/static/cmake/modules/FindNETCDF.cmake)

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
patch ./src/cmake/modules/FindGLIB.cmake < ./configs/${VER%.*}/static/cmake/modules/FindGLIB.cmake.patch
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
	-DBLAS_LIBRARIES="-lblas;${LIBGFORTRAN_ROOT}/libgfortran.a;${LIBGFORTRAN_ROOT}/libquadmath.a;${LIBGCC}" \
	-DCURL_INCLUDE_DIR="${CURL_ROOT}/include" \
	-DCURL_LIBRARY="${CURL_ROOT}/lib/libcurl.a" \
	-DGDAL_EXTRA_LIBS="${GDAL_EXTRA_LIBS}" \
	-DGLIB_INCLUDE_DIRS="${GLIB_H_DIR};${GLIBCONFIG_H_DIR}" \
	-DGLIB_LIBRARY="${LIBGLIB} ${LIBGLIB_ROOT}/libgthread-2.0.a ${LIBINTL} -liconv ${LIBPCRE} -framework Foundation" \
	-DLAPACK_LIBRARIES="-llapack;${LIBGFORTRAN_ROOT}/libgfortran.a;${LIBGFORTRAN_ROOT}/libquadmath.a;${LIBGCC}" \
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
