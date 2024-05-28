#!/bin/bash
set -eu

#Some cleanup
rm -rf install src
rm -rf Python-2.7.3
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh https://issm.ess.uci.edu/files/externalpackages/Python-2.7.3.tgz Python-2.7.3.tgz

#Untar and move python into install directory
tar -zxvf  Python-2.7.3.tgz
mv Python-2.7.3/* src
rm -rf Python-2.7.3

#Configure python
cd src 
./configure \
 --prefix="$ISSM_DIR/externalpackages/python/install" \
 --enable-shared

if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install

cd ../install/include
cp python2.7/* ./
cd ../lib
ln -s  libpython2.7.so.1.0 libpython.so
