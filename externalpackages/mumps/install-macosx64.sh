#!/bin/bash
set -eu

#Some cleanup
rm -rf src install mumps-4.10.0-p3
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/mumps-4.10.0-p3.tar.gz' 'mumps-4.10.0-p3.tar.gz'

#Untar 
tar -zxvf  mumps-4.10.0-p3.tar.gz

#Move mumps into src directory
mv mumps-4.10.0-p3/* src
rm -rf mumps-4.10.0-p3

#configuration: 
cp configs/Makefile-macosx64.inc src/Makefile.inc

#Configure mumps
cd src

#Compile mumps
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
