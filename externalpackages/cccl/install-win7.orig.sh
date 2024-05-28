#!/bin/bash
set -eu

#Some cleanup
rm -rf src install cccl-0.03
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/cccl-0.03.tar.gz' 'cccl-0.03.tar.gz'

#Untar 
tar -zxvf  cccl-0.03.tar.gz

#Move cccl into install directory
mv cccl-0.03/* src
rm -rf cccl-0.03

cd src 

#Compile
./configure --prefix="$ISSM_DIR/externalpackages/cccl/install"

make
make install
