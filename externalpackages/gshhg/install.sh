#!/bin/bash
set -eu


## Constants
#
VER="2.3.7"

PREFIX="${ISSM_DIR}/externalpackages/gshhg/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX}
mkdir -p ${PREFIX}

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-gmt-${VER}.tar.gz" "gshhg-gmt-${VER}.tar.gz"

# Unpack source
tar -zxvf gshhg-gmt-${VER}.tar.gz

# Install
mv gshhg-gmt-${VER}/* ${PREFIX}
