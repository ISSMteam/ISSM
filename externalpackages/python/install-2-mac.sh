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

#Configure and compile
cd src 
./configure \
 --enable-framework="$ISSM_DIR/externalpackages/python/install/Library/Frameworks" 
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install

cd ../install

#get rid of bin, because it's just a copy of
#Library/Frameworks/Python.framework/Versions/2.7/bin, and will not reflect
#new changes being made
rm -rf bin
ln -s Library/Frameworks/Python.framework/Headers include
ln -s Library/Frameworks/Python.framework/Versions/2.7/lib lib
ln -s Library/Frameworks/Python.framework/Versions/2.7/bin bin

#Patch pyport.h:
cd include
patch pyport.h $ISSM_DIR/externalpackages/python/patches/pyport.h.patch
