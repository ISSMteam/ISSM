#!/bin/bash
set -eu

#Some cleanup
rm -rf install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh https://issm.ess.uci.edu/files/externalpackages/modules-3.2.9c.tar.gz modules-3.2.9c.tar.gz

#Untar and move python into install directory
tar -zxvf  modules-3.2.9c.tar.gz
mv modules-3.2.9 install

#Configure doxygen
cd install
./configure \
  --prefix "$ISSM_DIR/externalpackages/modules/install" \
  --without-x

#compile and install
make
make install
