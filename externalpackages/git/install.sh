#!/bin/bash
set -eu

#Some cleanup
rm -rf src install git-1.7.10.2
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/git-1.7.10.2.tar.gz' 'git-1.7.10.2.tar.gz'

#Untar 
tar -zxvf  git-1.7.10.2.tar.gz

#Move git into install directory
mv git-1.7.10.2/* src
rm -rf git-1.7.10.2

#install
cd src 
./configure  --prefix="$ISSM_DIR/externalpackages/git/install" 
#--with-python="$ISSM_DIR/externalpackages/python/install/bin/python" #Do we really need this line?
	
#Compile
make install
