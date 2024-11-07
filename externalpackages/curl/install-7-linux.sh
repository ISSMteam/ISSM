#!/bin/bash
set -eu


## Constants
#
VER="7.73.0"

PREFIX="${ISSM_DIR}/externalpackages/curl/install" # Set to location where external package should be installed

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://curl.se/download/curl-${VER}.tar.gz" "curl-${VER}.tar.gz"

# Unpack source
tar -zxvf curl-${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Move source to 'src' directory
mv curl-${VER}/* src
rm -rf curl-${VER}

# Configure
cd src
./configure \
	--prefix="${PREFIX}" \
	--disable-static \
	--disable-dependency-tracking \
	--disable-manual \
	--disable-verbose \
	--with-zlib="${ZLIB_ROOT}" \
	--without-ssl

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
