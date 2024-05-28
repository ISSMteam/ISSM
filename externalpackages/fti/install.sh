#!/bin/bash
set -eu

#Some cleanup
rm -rf install fti-0.9.2
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/fti-0.9.2.tgz' 'fti-0.9.2.tgz'

#Untar 
tar -zxvf  fti-0.9.2.tgz

#Move mpich into src directory
mv fti/* install
rm -rf fti

#Configure mpich
cd install
cp ../Makefile .

#Compile mpich (this new version supports parallel make)
make
make install 
