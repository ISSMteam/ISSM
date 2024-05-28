#!/bin/bash
set -eu

#Some cleanup
rm -rf src install math77
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/math77.tar.gz' 'math77.tar.gz'

#Untar 
tar -zxvf  math77.tar.gz

#Move math77 into src directory
mv math77/* src
rm -rf math77

#Configure math77
cd src

#Compile math77
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
