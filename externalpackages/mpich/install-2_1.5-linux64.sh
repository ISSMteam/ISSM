#!/bin/bash
set -eu

#Some cleanup
rm -rf src install mpich2-1.5
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/mpich2-1.5.tar.gz' 'mpich2-1.5.tar.gz'

#Untar 
tar -zxvf  mpich2-1.5.tar.gz

#Move mpich2 into src directory
mv mpich2-1.5/* src
rm -rf mpich2-1.5

#Configure mpich2
cd src
./configure \
	--prefix="$ISSM_DIR/externalpackages/mpich/install" \
	--enable-shared \
	--enable-sharedlibs=gcc \
	--enable-f91=gfortran  CC=gcc FC=gfortran

#Compile mpich2 (this new version supports parallel make)
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
