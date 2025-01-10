#!/bin/bash
set -e


## TODO
#	- May want to supply path to Python instead of, effectively, using result of `which python`
#

## Constants
#
VER="3.10.0"

PREFIX="${ISSM_DIR}/externalpackages/gdal/install"

## Environment
#
export CFLAGS="${CFLAGS} -w"
export CXXFLAGS="${CXXFLAGS} -w"

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://github.com/OSGeo/gdal/releases/download/v${VER}/gdal-${VER}.tar.gz" "gdal-${VER}.tar.gz"

# Unpack source
tar -zxvf gdal-${VER}.tar.gz

# Move source into 'src' directory
mv gdal-${VER}/* src
rm -rf gdal-${VER}

# Configure
cd src
cmake \
	-DCMAKE_INSTALL_PREFIX="${PREFIX}" \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DBUILD_SHARED_LIBS=ON \
	-DBUILD_PYTHON_BINDINGS=ON \
	-DGDAL_PYTHON_INSTALL_PREFIX="${PREFIX}" \
	-DGDAL_SET_INSTALL_RELATIVE_RPATH=ON \
	-DGDAL_USE_ARCHIVE=OFF \
	-DGDAL_USE_CRYPTOPP=OFF \
	-DCURL_DIR="${CURL_ROOT}" \
	-DGDAL_USE_JPEG_INTERNAL=ON \
	-DGDAL_USE_JPEG12_INTERNAL=ON \
	-DNETCDF_DIR="${NETCDF_ROOT}" \
	-DGDAL_USE_OPENJPEG=OFF \
	-DGDAL_USE_OPENSSL=OFF \
	-DGDAL_USE_PNG_INTERNAL=ON \
	-DPROJ_DIR="${PROJ_ROOT}" \
	-DGDAL_USE_TIFF_INTERNAL=ON \
	-DZLIB_DIR="${ZLIB_ROOT}" \
	-DGDAL_USE_ZSTD=OFF

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
