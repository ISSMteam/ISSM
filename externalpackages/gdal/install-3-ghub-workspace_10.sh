#!/bin/bash
set -eu


## TODO
#	- May want to supply path to Python instead of, effectively, using result of `which python`
#

## Constants
#
VER="3.10.0"

## Environment
#
export PREFIX="${ISSM_EXT_DIR}/gdal/install" # NOTE: Need to export this to properly set destination root for Python libraries on macOS (should not affect Linux build). Set to location where external package should be installed.

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
./configure \
	--prefix="${PREFIX}" \
	--enable-fast-install \
	--with-python="/apps/share64/debian10/anaconda/anaconda-7/bin/python3.8" \
	--with-curl="${CURL_ROOT}/bin/curl-config" \
	--with-hdf5="${HDF5_ROOT}" \
	--with-libz="${ZLIB_ROOT}" \
	--with-netcdf="${NETCDF_ROOT}" \
	--with-proj="${PROJ_ROOT}"

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
