#!/bin/bash
set -eu 

#Some cleanup
rm -rf src
rm -rf ipython-1.0.0
mkdir src 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/ipython-1.0.0.tar.gz' 'ipython-1.0.0.tar.gz'

#Untar 
tar -zxvf  ipython-1.0.0.tar.gz

#Move ipython into src directory
mv ipython-1.0.0/* src
rm -rf ipython-1.0.0

#install  ipython
cd src
python setup.py build
python setup.py install
