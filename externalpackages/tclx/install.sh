#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf tclx8.4
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/tclx8.4.tar.gz' 'tclx8.4.tar.gz'

#Untar 
tar -zxvf  tclx8.4.tar.gz

#Move tclx into src directory
mv tclx8.4/* src
rm -rf tclx8.4

#Configure tclx
cd src
./configure --prefix="$ISSM_DIR/externalpackages/tclx/install"  \
			--exec-prefix="$ISSM_DIR/externalpackages/tclx/install"  \
	        --with-tcl=$ISSM_DIR/externalpackages/tcl/install/lib

#Compile and install tclx
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
