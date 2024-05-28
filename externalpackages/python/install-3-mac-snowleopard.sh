#!/bin/bash
set -eu

#Some cleanup
rm -rf install src
rm -rf Python-3.2.2
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh https://issm.ess.uci.edu/files/externalpackages/Python-3.2.2.tgz Python-3.2.2.tgz

#exports
export CC
export MACOSX_DEPLOYMENT_TARGET=10.6

#Untar and move python into install directory
tar -zxvf  Python-3.2.2.tgz
mv Python-3.2.2/* src
rm -rf Python-3.2.2

#Configure doxygen
cd src 
# --enable-framework needs to have the form "$SOME_PATH/Library/Frameworks" to avoid installing components in /Applications directory
# --prefix is recognized as $SOME_PATH as long as this form is taken, so it's not necessary to include
./configure \
 --enable-framework="$ISSM_DIR/externalpackages/python/install/Library/Frameworks" 

#make
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
cd ..

cd install/bin
ln -s python3.2 python 
cd ../
ln -s Python.framework/Versions/3.2/include/python3.2m include
ln -s Python.framework/Versions/3.2/lib/ lib
