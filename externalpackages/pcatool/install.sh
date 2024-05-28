#!/bin/bash

#Some cleanup
rm -rf install pcatool
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/pcatool.tar.gz' 'pcatool.tar.gz'

#Untar  into install
cd install 
tar -zxvf  ../pcatool.tar.gz
