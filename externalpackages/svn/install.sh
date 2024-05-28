#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf subversion-1.6.18
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/subversion-1.6.18.tar.gz' 'subversion-1.6.18.tar.gz'

#Untar 
tar -zxvf  subversion-1.6.18.tar.gz

#Move subversion into src directory
mv subversion-1.6.18/* src
rm -rf subversion-1.6.18

#Configure subversion
cd src
./configure --prefix="$ISSM_DIR/externalpackages/svn/install" \
	--with-swig="$ISSM_DIR/externalpackages/swig/install"  \
	PYTHON2="$ISSM_DIR/externalpackages/python/install/bin/python" \
	--with-sqlite="$ISSM_DIR/externalpackages/sqlite/install"

#Compile and install subversion
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
