#!/bin/bash
set -eu

#Erase install
rm -rf install  src ColPack

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/ColPack-1.0.3.tar.gz' 'ColPack-1.0.3.tar.gz'

#install directory
mkdir src
tar -zxvf ColPack-1.0.3.tar.gz 
mv ColPack/* src
rm -rf ColPack

#compile
cd src
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
cd ..

#install
ln -s src/build ./install
