#!/bin/bash
set -eu

#Some cleanup
rm -rf src install lapack-3.4.1 lapack-3.4.1.tgz
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/lapack-3.4.1.tgz' 'lapack-3.4.1.tgz'

#Untar 
tar -zxvf  lapack-3.4.1.tgz

#Move lapack into install directory
mv lapack-3.4.1/* src
rm -rf lapack-3.4.1

#install
cd src 
cp ../configs/macosx64/make.inc ./

#Compile and install lapack
if [ $# -eq 0 ]; then
	make lib
else
	make -j $1 lib
fi

#Compile 
cd ../install
mkdir lib
cd lib
cp ../../src/liblapack.a .
