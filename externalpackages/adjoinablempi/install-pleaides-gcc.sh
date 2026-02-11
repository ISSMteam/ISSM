#!/bin/bash
set -eu


# Cleanup
rm -rf install src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://github.com/ISSMteam/ExternalPackages/raw/refs/heads/main/adjoinablempi.tar.gz' 'adjoinablempi.tar.gz'

# Unpack source
tar -zxf adjoinablempi.tar.gz

# Configure
cd src
./configure \
	--prefix="$ISSM_DIR/externalpackages/adjoinablempi/install" \
	--libdir="$ISSM_DIR/externalpackages/adjoinablempi/install/lib" \
	--with-mpi-root="$ISSM_DIR/externalpackages/mpich/install" \
	--enable-requestOnTrace CFLAGS="-g -O0"

# Compile
make clean
if [ $# -eq 0 ]; then
	make 
else
	make -j $1
fi
make install
