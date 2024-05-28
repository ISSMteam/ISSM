#!/bin/bash
set -eu

#Some cleanup
rm -rf install m2html
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/m2html.zip' 'm2html.zip'

#uncompress
unzip m2html.zip

#move to install directory
mv m2html/* install
rm -rf m2html

#patch m2html
cd install
patch m2html.m ../m2html.m.patch
cd ..
