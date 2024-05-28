#!/bin/bash
set -eu


## Environment
#
export LD_LIBRARY_PATH="" # Ensure that libtool does not hardcode local paths set by running $ISSM_DIR/etc/environment.sh into binaries
export LD_RUN_PATH="" # Ensure that libtool does not hardcode local paths set by running $ISSM_DIR/etc/environment.sh into binaries

## Constants
#
VER="3.2"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/mpich-${VER}.tar.gz" "mpich-${VER}.tar.gz"

# Unpack source
tar -zxvf mpich-${VER}.tar.gz

# Cleanup
rm -rf install src
mkdir install src

# Move source into 'src' directory
mv mpich-${VER}/* src
rm -rf mpich-${VER}

# Configure
cd src
./configure \
	--prefix="${ISSM_DIR}/externalpackages/mpich/install" \
	--disable-shared \
	--disable-dependency-tracking \
	--enable-fast=all \
	--with-pic

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi

# Return to initial directory
cd ..

# Strip RPATH/RUNPATH entries from executables
#
# NOTE:
# - We are doing this so that we do not ship executables that have hardcoded
#	local paths in their RPATH/RUNPATH entries
# - This might otherwise be accomplished with extensive changes to libtool's
#	handling of rpath
#
chrpath -d ./install/bin/hydra_pmi_proxy
chrpath -d ./install/bin/mpiexec
