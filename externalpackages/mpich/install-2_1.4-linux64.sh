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
export -n F90 
export CFLAGS="$CFLAGS -fPIC"
export FFLAGS="$FFLAGS -fPIC"
./configure \
	--prefix="$ISSM_DIR/externalpackages/mpich/install" \
	--enable-shared \
	--enable-sharedlibs=gcc \
	--enable-f91=gfortran 

#Compile mpich2 (parallel make not supported)
make
make install 
