#!/bin/bash
set -eu


## Constants
#
VER="1.0.0"

PREFIX="${ISSM_DIR}/externalpackages/shell2junit/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX}
mkdir -p ${PREFIX}/bin

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/shell2junit-${VER}.zip" "shell2junit-${VER}.zip"

# Unpack source
unzip -q shell2junit-${VER}.zip

# Install
mv shell2junit-${VER}/sh2ju.sh ${PREFIX}/bin

# Cleanup
rm -rf shell2junit-${VER}
