#!/bin/bash
set -eu

#Some cleanup
rm -rf install src make-3.82
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/make-3.82.tar.gz' 'make-3.82.tar.gz'

#Untar 
tar -zxvf  make-3.82.tar.gz

#Move make into install directory
mv make-3.82/* src
rm -rf make-3.82

#Apply patches
cd src 

#Configure and compile: 
./configure --prefix=$ISSM_DIR/externalpackages/gmake/install

if [ $# -eq 0 ]; then
	make
	make install
else 
	make -j $1; 
	make -j $1 install;
fi
