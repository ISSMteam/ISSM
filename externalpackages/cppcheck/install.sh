#!/bin/bash
set -eu

#Erase install
rm -rf install  src cppcheck-1.48
mkdir src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/cppcheck-1.48.tar' 'cppcheck-1.48.tar'
tar -xvf cppcheck-1.48.tar

mv cppcheck-1.48/* src
rm -rf cppcheck-1.48

#compile
cd src
if [ $# -eq 0 ]; then
	make 
else 
	make -j $1
fi  
make install PREFIX=$ISSM_DIR/externalpackages/cppcheck/install
cd ..
