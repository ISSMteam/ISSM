#!/bin/bash
set -eu

#Some cleanup
rm -rf install
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/tex2im-1.8.tar.gz' 'tex2im-1.8.tar.gz'

#Untar 
tar -zxvf  tex2im-1.8.tar.gz

#Move tex2im into src directory
mv tex2im-1.8/tex2im install/
rm -rf tex2im-1.8
