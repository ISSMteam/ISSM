#!/bin/bash
set -eu

#Some cleanup
rm -rf install
rm -rf modules-3.2.9c
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/modules-3.2.9c.tar.gz' 'modules-3.2.9c.tar.gz'

#Untar 
tar -zxvf  modules-3.2.9c.tar.gz

#Move modules into src directory

#Configure modules
cd modules-3.2.9
./configure \
	--prefix=$ISSM_DIR/externalpackages/modules/install\
	--with-tcl-lib=$ISSM_DIR/externalpackages/tcl/install/lib\
	--with-tcl-inc=$ISSM_DIR/externalpackages/tcl/install/include\
	--with-tcl-ver=8.5 \
	--with-tclx-lib=$ISSM_DIR/externalpackages/tclx/install/lib/tclx8.4\
    --with-tclx-inc=$ISSM_DIR/externalpackages/tclx/install/include\
	--with-tclx-ver=8.4 \
	--with-version-path=/usr/local/modules/versions \
	--with-skel-path=/usr/local/modules/etc/skel \
	--with-etc-path=/usr/local/modules/etc \
	--with-module-path=/usr/local/modules/files \
	--disable-dependency-tracking

#Compile and install modules
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
sudo make install
