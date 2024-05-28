#!/bin/bash
set -eu


## Constants
#
VER="7.73.0"

PREFIX="${ISSM_DIR}/externalpackages/curl/install" # Set to location where external package should be installed

## Environment
#
export MACOSX_DEPLOYMENT_TARGET="10.5" # Allows fall back to older API (source: https://curl.se/docs/install.html)

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/curl-${VER}.tar.gz" "curl-${VER}.tar.gz"

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
	--disable-shared \
	--disable-dependency-tracking \
	--disable-manual \
	--disable-verbose \
	--disable-ldap \
	--disable-ldaps \
	--with-zlib="${ZLIB_ROOT}" \
	--without-zstd \
	--without-libidn2 \
	--without-nghttp2 \
	--without-brotli \
	--without-librtmp \
	--without-ssl

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
