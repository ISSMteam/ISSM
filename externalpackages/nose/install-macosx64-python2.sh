#!/bin/bash
#Install Python nose module

rm -rf src  install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/nose-1.1.2.tar.gz' 'nose-1.1.2.tar.gz'
tar -zxvf  nose-1.1.2.tar.gz
mv nose-1.1.2 src
rm -rf nose-1.1.2

cd src
python ./setup.py build
python ./setup.py install

#to be flagged by jenkins, we create an empty install dir: 
cd ../
mkdir install
touch install/emptyfile
