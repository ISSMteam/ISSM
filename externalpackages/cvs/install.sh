#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf cvs-1.11.23
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/cvs-1.11.23.tar.gz' 'cvs-1.11.23.tar.gz'

#Untar 
tar -zxvf  cvs-1.11.23.tar.gz

#Move subversion into src directory
mv cvs-1.11.23/* src
rm -rf cvs-1.11.23

#Patch getline
cd src
cat lib/getline.c | sed -e "s/getline /get_line /" > BACK && mv BACK lib/getline.c
cat lib/getline.h | sed -e "s/getline /get_line /" > BACK && mv BACK lib/getline.h

#Configure
./configure --prefix="$ISSM_DIR/externalpackages/cvs/install" 

#Compile and install subversion
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
