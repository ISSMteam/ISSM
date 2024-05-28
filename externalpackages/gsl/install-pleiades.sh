#!/bin/bash
set -eu

#Some cleanup
rm -rf src install gsl-1.15
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/gsl-1.15.tar.gz' 'gsl-1.15.tar.gz'

#Untar 
tar -zxvf  gsl-1.15.tar.gz

#Move gsl into src directory
mv gsl-1.15/* src
rm -rf gsl-1.15

#Configure gsl
cd src

./configure \
	--prefix="$ISSM_DIR/externalpackages/gsl/install" 

#Compile gsl
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
