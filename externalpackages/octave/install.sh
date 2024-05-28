#!/bin/bash
set -eu

#Some cleanup
rm -rf src install octave-3.6.2 
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/octave-3.6.2.tar.gz" "octave-3.6.2.tar.gz"

#Untar 
tar -zxvf  octave-3.6.2.tar.gz

#Move octave into install directory
mv octave-3.6.2/* src
rm -rf octave-3.6.2

#install
cd src 
./configure \
 --prefix=$ISSM_DIR/externalpackages/octave/install \
 --disable-readline

if [ $# -eq 0 ];
then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
