#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf pysvn-1.7.6
mkdir src 

export CC=/usr/local/gfortran/bin/gcc
export CXX=/usr/local/gfortran/bin/g++
export FFLAGS=-ff2c

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/pysvn-1.7.6.tar.gz' 'pysvn-1.7.6.tar.gz'

#Untar 
tar -zxvf  pysvn-1.7.6.tar.gz

#Move pysvn into src directory
mv pysvn-1.7.6/* src
rm -rf pysvn-1.7.6

$install
cd src
python setup.py build
python setup.py install
