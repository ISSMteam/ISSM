#!/bin/bash

#Some cleanup
rm -rf src install mpich2-1.4
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/mpich2-1.4.tar.gz' 'mpich2-1.4.tar.gz'

#Untar 
tar -zxvf  mpich2-1.4.tar.gz

#Move mpich2 into src directory
mv mpich2-1.4/* src
rm -rf mpich2-1.4

#Configure mpich2
cd src
export FCFLAGS=" -m64"
export FFLAGS=" -m64"
export CFLAGS=" -arch x86_64"
export CXXFLAGS=" -arch x86_64"
./configure \
	--prefix="$ISSM_DIR/externalpackages/mpich/install" \
	--enable-f91 \
	--enable-sharedlibs=osx-gcc \
	--enable-shared \
	--enable-fc

#Compile mpich2 (parallel make not supported)
make
make install 
