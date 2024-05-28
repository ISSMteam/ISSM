#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf apr-util-1.4.1
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/apr-util-1.4.1.tar.gz' 'apr-util-1.4.1.tar.gz'

#Untar 
tar -zxvf  apr-util-1.4.1.tar.gz

#Move apr-util into src directory
mv apr-util-1.4.1/* src
rm -rf apr-util-1.4.1

#Configure apr-util
cd src
./configure  --prefix="$ISSM_DIR/externalpackages/apr-util/install" --with-apr="$ISSM_DIR/externalpackages/apr/install"

#Compile and install apr-util
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
