#!/bin/bash
set -eu

#Some cleanup
rm -rf install
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/ncview-1.93g.tar.gz' 'ncview-1.93g.tar.gz'

#Untar 
tar -zxvf  ncview-1.93g.tar.gz

#Move doxygen into src directory
mv ncview-1.93g/* install
rmdir ncview-1.93g

#Configure doxygen
cd install
./configure \
	--prefix="$ISSM_DIR/externalpackages/ncview/install" \
	--x-libraries=/usr/X11/lib/ \
	--x-includes=/usr/X11/include/ \
	--with-netcdf_incdir="$ISSM_DIR/externalpackages/netcdf/install/include/" \
	--with-netcdf_libdir="$ISSM_DIR/externalpackages/netcdf/install/lib/"

make
make install
