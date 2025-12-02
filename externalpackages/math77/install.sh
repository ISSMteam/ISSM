#!/bin/bash
set -eu


## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/math77/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX}

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://github.com/ISSMteam/ExternalPackages/raw/refs/heads/main/math77.tar.gz" "math77.tar.gz"

# Unpack source
tar -zxvf math77.tar.gz

# Move math77/ to/ src/ directory
mv math77 src

# Compile and install
cd src
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
