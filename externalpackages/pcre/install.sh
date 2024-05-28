#!/bin/bash
set -eu

#Cleaning
rm -rf install
rm -rf pcre-8.21
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/pcre-8.21.tar.gz' 'pcre-8.21.tar.gz'

#Untar and move python into install directory
tar -zxvf  pcre-8.21.tar.gz
mv pcre-8.21/* install
rm -rf pcre-8.21

#Configure doxygen
cd install 
./configure --prefix "$ISSM_DIR/externalpackages/python/install"
make
make install
