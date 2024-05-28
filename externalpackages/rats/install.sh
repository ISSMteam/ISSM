#!/bin/bash
set -eu

#Erase install
rm -rf install  src rats-2.3

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/rats-2.3.tar.gz' 'rats-2.3.tar.gz'

#install directory
mkdir src
tar -zxvf rats-2.3.tar.gz 
mv rats-2.3/* src
rm -rf rats-2.3

#compile
cd src
./configure --prefix=$ISSM_DIR/externalpackages/rats/install
make
make install
cd ..
