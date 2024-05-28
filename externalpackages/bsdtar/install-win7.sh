#!/bin/bash
set -eu

#Some cleanup
rm -rf install src libarchive-3.0.3
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/libarchive-3.0.3.tar.gz' 'libarchive-3.0.3.tar.gz'

#Untar 
tar -zxvf  libarchive-3.0.3.tar.gz

#Move libarchive into src directory
mv libarchive-3.0.3/* src
rm -rf libarchive-3.0.3

cd src 
./configure --prefix="$ISSM_DIR/externalpackages/bsdtar/install" 
make
make install
