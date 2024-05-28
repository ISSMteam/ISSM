#!/bin/bash
set -eu

#Some cleanup
rm -rf install xerces-c-src_2_8_0 src
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/xerces-c-src_2_8_0.tar.gz' 'xerces-c-src_2_8_0.tar.gz'

#Untar 
tar -zxvf  xerces-c-src_2_8_0.tar.gz

#Move xerces-c-tools into install directory
mv xerces-c-src_2_8_0/* src
rm -rf xerces-c-src_2_8_0

#Apply patches
cd src/src/xercesc/

#Configure
./runConfigure -plinux -cgcc -xg++ -minmem -nsocket -tnative -rnone -s 

#Compile xerces-c-tools
make
