#!/bin/bash
set -eu


## Constants
#
VER="2.3.4"

PREFIX="${ISSM_DIR}/externalpackages/gshhg/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX}
mkdir -p ${PREFIX}

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/gshhg-gmt-${VER}.tar.gz" "gshhg-gmt-${VER}.tar.gz"

# Unpack source
tar -zxvf gshhg-gmt-${VER}.tar.gz

# Install
mv gshhg-gmt-${VER}/* ${PREFIX}
