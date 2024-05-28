#!/bin/bash
set -eu

#Some cleanup
rm -rf install
rm -rf libermate-0.4
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/libermate-0.4.tgz' 'libermate-0.4.tgz'

#Untar 
tar -zxvf  libermate-0.4.tgz

#Move libermate into install directory
mv libermate-0.4/* install
rm -rf libermate-0.4
