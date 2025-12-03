#!/bin/bash
set -eu


## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/vim/install" # Set to location where external package should be installed

VER="9.1.1943"

# Cleanup
rm -rf ${PREFIX} src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://github.com/vim/vim/archive/refs/tags/v${VER}.tar.gz" "v${VER}.tar.gz"

# Unpack source
tar -zxvf v${VER}.tar.gz
mv vim-${VER} src

# Configure (icc seems to have issues with wctype.h)
export CC=gcc
cd src/src
./configure \
	--prefix="${PREFIX}" \
	--with-gcc="/usr/bin/gcc" \
	--with-tlib="/lib/"

# Compile and install
make
make install
