#!/bin/bash
set -eu

#Some cleanup
rm -rf install netcdf-cxx-4.2
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/netcdf-cxx-4.2.tar.gz' 'netcdf-cxx-4.2.tar.gz'

#Untar 
tar -zxvf netcdf-cxx-4.2.tar.gz

#Move metis into install directory
mv netcdf-cxx-4.2/* install
rm -rf netcdf-cxx-4.2

#Compile
export CXXFLAGS="-I$ISSM_DIR/externalpackages/netcdf/install/include "
cd install 
./configure --prefix="$ISSM_DIR/externalpackages/netcdf-cxx/install" 
make
make install
