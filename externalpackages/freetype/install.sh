#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf freetype-2.5.0
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/freetype-2.5.0.tar.gz' 'freetype-2.5.0.tar.gz'

#Untar 
tar -zxvf  freetype-2.5.0.tar.gz

#Move freetype into src directory
mv freetype-2.5.0/* src
rm -rf freetype-2.5.0

#Configure freetype
cd src
sudo ./configure 

#Compile and install freetype
if [ $# -eq 0 ]; then
	sudo make
else
	sudo make -j $1
fi
sudo make install
