#!/bin/bash
set -eu
#metis 5.0 should be used: srand48 and drand48 are being redefined in conflict to the stdlib equivalent functions.

#Some cleanup
rm -rf install metis-5.0.1
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/metis-5.0.1.tar.gz' 'metis-5.0.1.tar.gz'

#Untar 
tar -zxvf  metis-5.0.1.tar.gz

#Move metis into install directory
mv metis-5.0.1/* install
rm -rf metis-5.0.1

#Apply patches
cd install 

#Compile metis
make config prefix="$ISSM_DIR/externalpackages/metis/install"
make install
