#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf apr-1.4.6
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/apr-1.4.6.tar.gz' 'apr-1.4.6.tar.gz'

#Untar 
tar -zxvf  apr-1.4.6.tar.gz

#Move apr into src directory
mv apr-1.4.6/* src
rm -rf apr-1.4.6

#Configure apr
cd src
./configure  --prefix="$ISSM_DIR/externalpackages/apr/install" 

#Compile and install apr
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
