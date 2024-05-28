#!/bin/bash
set -eu


# NOTE:
# - This installation script should work equally well under Linux. If your 
#	Linux distribution does not include libtiff, please,
#		- copy this file to install-4-linux.sh
#		- verify that it successfully compiles and installs libtiff
#		- commit the new script (including any modifications) to the repo
#

# Constants
#
VER="4.2.0"

PREFIX="${ISSM_EXT_SHARED_DIR}/libtiff/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/libtiff-v${VER}.tar.gz" "libtiff-v${VER}.tar.gz"

# Unpack source
tar -zxvf libtiff-v${VER}.tar.gz

# Move source into 'src' directory
mv libtiff-v${VER}/* src
rm -rf libtiff-v${VER}

# Configure
cd src
cmake -S . -B ${PREFIX}
cd ${PREFIX}

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
