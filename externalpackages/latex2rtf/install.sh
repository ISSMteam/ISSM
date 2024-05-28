#!/bin/bash
set -eu

#Some cleanup
rm -rf latex2rtf-2.0.0 cfg install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/latex2rtf-2.0.0.tar.gz' 'latex2rtf-2.0.0.tar.gz'

#untar 
tar -zxvf  latex2rtf-2.0.0.tar.gz
mv latex2rtf-2.0.0 install

#Compile
cd install
export PREFIX="$ISSM_DIR/externalpackages/latex2rtf/install/"
make
