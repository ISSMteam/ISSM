#!/bin/bash
set -eu

#Some cleanup
rm -rf src install bbftp-client-3.2.0

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/bbftp-client-3.2.0.tar.gz' 'bbftp-client-3.2.0.tar.gz'

#Create install directories
mkdir install src

#Untar 
tar -zxvf  bbftp-client-3.2.0.tar.gz

#Move bbftp-client into install directory
mv bbftp-client-3.2.0/* src
rm -rf bbftp-client-3.2.0

#Apply patches
cd src 

#Configure and compile
cd bbftpc
./configure --prefix=$ISSM_DIR/externalpackages/bbftp/install
make
make install
