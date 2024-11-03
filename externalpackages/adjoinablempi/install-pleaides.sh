#!/bin/bash
set -eu


# Cleanup
rm -rf install src

# Download source
#$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/adjoinablempi.tar.gz' 'adjoinablempi.tar.gz'

# Unpack source
tar -zxf adjoinablempi.tar.gz

# Configure
cd src
./configure \
	--prefix="$ISSM_DIR/externalpackages/adjoinablempi/install" \
	--libdir="$ISSM_DIR/externalpackages/adjoinablempi/install/lib" \
	--with-mpi-root="/nasa/sgi/mpt/2.06rp16/" \
	--enable-requestOnTrace

# Compile
make clean
if [ $# -eq 0 ]; then
	make 
else
	make -j $1
fi
make install
