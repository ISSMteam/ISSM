#!/bin/bash
set -e


## Constants
#
VER="1.14.3"

PREFIX="${ISSM_DIR}/externalpackages/hdf5/install" # Set to location where external package should be installed

## Environment
#
export CC=mpicc
export CFLAGS="${CFLAGS} -w"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_14_3/src/hdf5-${VER}.tar.gz" "hdf5-${VER}.tar.gz"

# Untar source
tar -zxvf hdf5-${VER}.tar.gz

# Cleanup
rm -rf install src
mkdir install src

# Move source to 'src' directory
mv hdf5-${VER}/* src/
rm -rf hdf5-${VER}

# Configure
cd src
./configure \
	--prefix="${PREFIX}" \
	--disable-dependency-tracking \
	--disable-static \
	--with-zlib="${ZLIB_ROOT}" \
	--with-szlib="no" \
	--enable-hl

# Compile and install
#
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
