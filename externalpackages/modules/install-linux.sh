#!/bin/bash
set -eu


## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/modules/install" # Set to location where external package should be installed

VER="5.6.1"

# Cleanup
rm -rf ${PREFIX} src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://github.com/envmodules/modules/archive/refs/tags/v${VER}.tar.gz" "v${VER}.tar.gz"

# Unpack source
tar -zxvf v${VER}.tar.gz
mv modules-${VER} src

# Configure
cd src
./configure \
	--prefix "${PREFIX}" \
	--without-x

# Compile and install
make
make install
