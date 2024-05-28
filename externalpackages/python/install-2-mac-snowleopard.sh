#!/bin/bash
set -eu

#Some cleanup
rm -rf install src
rm -rf Python-2.7.3
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh https://issm.ess.uci.edu/files/externalpackages/Python-2.7.3.tgz Python-2.7.3.tgz

#exports
export MACOSX_DEPLOYMENT_TARGET=10.6

#Untar and move python into install directory
tar -zxvf  Python-2.7.3.tgz
mv Python-2.7.3/* src
rm -rf Python-2.7.3

#Configure doxygen
cd src 
# --enable-framework needs to have the form "$SOME_PATH/Library/Frameworks" to avoid installing components in /Applications directory
# --prefix is recognized as $SOME_PATH as long as this form is taken, so it's not necessary to include
./configure --enable-framework="$ISSM_DIR/externalpackages/python/install/Library/Frameworks" \
	LDFLAGS="-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin10/4.6.2/ -lgfortran "

#make
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install

cd ../install
ln -s Library/Frameworks/Python.framework/Headers include
ln -s Library/Frameworks/Python.framework/Versions/2.7/lib lib

#Patch pyport.h:
cd include
patch pyport.h $ISSM_DIR/externalpackages/python/patches/pyport.h.patch
