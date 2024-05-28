#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf libpng-1.5.10
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/libpng-1.5.10.tar.gz' 'libpng-1.5.10.tar.gz'

#Untar 
tar -zxvf  libpng-1.5.10.tar.gz

#Move libpng into src directory
mv libpng-1.5.10/* src
rm -rf libpng-1.5.10

#Configure libpng
cd src
sudo ./configure 

#Compile and install libpng
if [ $# -eq 0 ]; then
	sudo make
else
	sudo make -j $1
fi
sudo make install
