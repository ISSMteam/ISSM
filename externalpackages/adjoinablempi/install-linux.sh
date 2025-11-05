#!/bin/bash
set -eu


# Cleanup
rm -rf install src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://github.com/ISSMteam/ExternalPackages/raw/refs/heads/main/adjoinablempi.tar.gz' 'adjoinablempi.tar.gz'

# Unpack source
tar -zxvf adjoinablempi.tar.gz

# Configure
cd src
./configure \
	--prefix="${ISSM_DIR}/externalpackages/adjoinablempi/install" \
	--libdir="${ISSM_DIR}/externalpackages/adjoinablempi/install/lib" \
	--with-mpi-root="${ISSM_DIR}/externalpackages/petsc/install" \
	--enable-requestOnTrace

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
