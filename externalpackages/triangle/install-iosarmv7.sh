#!/bin/bash
set -eu

host_triplet="armv7-apple-darwin14.0.0"

#cleanup
rm -rf install triangle
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/triangle.zip' 'triangle.zip'

#untar
cd install
cp ../triangle.zip ./
unzip triangle.zip
rm triangle.zip

#copy new makefile
cp ../configs/ios/configure.make ./
cp ../makefile ./

#compile triangle
make

#patch triangle.h
patch triangle.h ../triangle.h.patch
