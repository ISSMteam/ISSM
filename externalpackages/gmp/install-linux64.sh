#!/bin/bash
set -eu

#Some cleanup
rm -rf src install gmp-5.0.5 
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/gmp-5.0.5.tar.bz2' 'gmp-5.0.5.tar.bz2'

#Untar 
bunzip2 gmp-5.0.5.tar.bz2
tar -xvf  gmp-5.0.5.tar

#Move gmp into install directory
mv gmp-5.0.5/* src
rm -rf gmp-5.0.5

#install
cd src 
./configure --prefix="$ISSM_DIR/externalpackages/gmp/install" 

#Compile and install gmp
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
