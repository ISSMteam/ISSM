#!/bin/bash
set -eu


## Environment
#
#export CPATH=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/malloc:/usr/include # If path to malloc is not included in shell configuration files, export it here after verifying that the path is correct.
export CFLAGS=-Wno-error=implicit-function-declaration

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
	--disable-dependency-tracking \
	--with-mpi-root="${ISSM_DIR}/externalpackages/petsc/install" \
	--enable-requestOnTrace

# Clean
make clean

# Compile
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi

# Install
make install
