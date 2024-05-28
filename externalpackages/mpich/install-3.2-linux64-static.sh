#!/bin/bash
set -eu

#Some cleanup
rm -rf src install mpich-3.2
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/mpich-3.2.tar.gz' 'mpich-3.2.tar.gz'

#Untar
tar -zxvf  mpich-3.2.tar.gz

#Move mpich into src directory
mv mpich-3.2/* src
rm -rf mpich-3.2

#Configure mpich
cd src
./configure \
	--prefix="$ISSM_DIR/externalpackages/mpich/install" \
	--disable-shared \
	--enable-strict=all \
	--enable-fast \
	--with-pic

#Compile mpich (this new version supports parallel make)
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi

make install
