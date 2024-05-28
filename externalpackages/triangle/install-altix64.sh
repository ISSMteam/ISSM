#!/bin/bash
set -eu

#Some cleanup 
rm -rf install triangle
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/triangle.zip' 'triangle.zip'

#Untar 
cd install
cp ../triangle.zip ./
unzip triangle.zip

#copy new makefile
cp ../configs/altix64/configure.make ./
cp ../makefile ./

#Compile triangle
make

#Patch triangle.h
patch triangle.h ../triangle.h.patch
