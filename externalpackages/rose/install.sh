#!/bin/bash
set -eu

#Some cleanup
rm -rf source build install
mkdir install source build

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/rose-0.9.5a-13219.tar.gz' 'rose-0.9.5a-13219.tar.gz'

#Untar 
tar -zxvf  rose-0.9.5a-13219.tar.gz

#Move rose into install directory
mv rose-0.9.5a-13219/* source
rm -rf rose-0.9.5a-13219

#Configure
cd build
../source/configure \
	--prefix=$ISSM_DIR/externalpackages/rose/install \
	--with-boost=$ISSM_DIR/externalpackages/boost/install\
	--srcdir=$ISSM_DIR/externalpackages/rose/source

if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
