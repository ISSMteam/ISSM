#!/bin/bash
set -eu

#Some cleanup
rm -rf install
rm -rf pyclips-1.0.7.348
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/pyclips-1.0.7.348.tar.gz' 'pyclips-1.0.7.348.tar.gz'

#Untar 
tar -zxvf  pyclips-1.0.7.348.tar.gz

#Move pyclips into install directory
mv pyclips/* install
rm -rf pyclips

#install
cd install
python setup.py build
sudo python setup.py install
