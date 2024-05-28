#!/bin/bash
set -eu

#Some cleanup
rm -rf install metis-4.0
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/metis-4.0.tar.gz' 'metis-4.0.tar.gz'

#Untar 
tar -zxvf  metis-4.0.tar.gz

#Move metis into install directory
mv metis-4.0/* install
rm -rf metis-4.0

#Apply patches
cd install 
patch -p1 < ../metis-4.0.patch
patch Makefile.in ../configs/4.0/cosmos/Makefile.in.patch

#Compile
make
