#!/bin/bash
set -eu

#Some cleanup
rm -rf install src

#Download latest version
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/openssl-0.9.8x.tar.gz' 'openssl-0.9.8x.tar.gz'

#Untar
tar -xzf openssl-0.9.8x.tar.gz
mv openssl-0.9.8x src
mkdir install

#Configure openssl
cd src
./config --prefix="$ISSM_DIR/externalpackages/openssl/install" shared
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
