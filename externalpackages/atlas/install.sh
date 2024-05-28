#!/bin/bash
set -eu

#Some cleanup
rm -rf src install ATLAS build
mkdir src install build

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/atlas3.10.0.tar.bz2' 'atlas3.10.0.tar.bz2'

#Untar 
tar -zxvf  atlas3.10.0.tar.bz2

#Move atlas into src directory
mv ATLAS/* src
rm -rf ATLAS

#Configure atlas
cd build

export CFLAGS=" -arch x86_64"

../src/configure \
	--prefix="$ISSM_DIR/externalpackages/atlas/install" 

#Compile atlas
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
