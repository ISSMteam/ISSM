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
	--with-tcl-lib="${ISSM_DIR}/externalpackages/tcl/install/lib" \
	--with-tcl-inc="${ISSM_DIR}/externalpackages/tcl/install/include" \
	--with-tcl-ver=8.5 \
	--with-tclx-lib="${ISSM_DIR}/externalpackages/tclx/install/lib/tclx8.4" \
	--with-tclx-inc="${ISSM_DIR}/externalpackagestclx/install/include" \
	--with-tclx-ver=8.4 \
	--with-version-path=/usr/local/modules/versions \
	--with-skel-path=/usr/local/modules/etc/skel \
	--with-etc-path=/usr/local/modules/etc \
	--with-module-path=/usr/local/modules/files \
	--disable-dependency-tracking

# Compile and install
make
sudo make install
