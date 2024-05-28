#!/bin/bash
set -eu

#Some cleanup
rm -rf src install qhull-2003.1
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/qhull-2003.1.tar.gz' 'qhull-2003.1.tar.gz'

#Untar 
tar -zxvf  qhull-2003.1.tar.gz

#Move qhull to src directory
rm -rf src/*
mv qhull-2003.1/* src/
rm -rf qhull-2003.1

#Configure qhull
cd src
./configure --prefix="$ISSM_DIR/externalpackages/qhull/install"
make
make install
