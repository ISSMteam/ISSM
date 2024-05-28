#!/bin/bash
set -eu

#Note of caution:  stop after boostrap phase, and run 
#bjam --debug-configuration, to figure out which paths boost is using to include 
#python. make sure everyone of these paths is covered by python. If not, just make 
#symlinks in externalpackages/python to what boost is expecting. Ther is NO WAY 
#to get the boost library to include python support without doing that. 

#Some cleanup
rm -rf install boost_1_55_0 src
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/boost_1_55_0.tar.gz' 'boost_1_55_0.tar.gz'

#Untar 
tar -zxvf  boost_1_55_0.tar.gz

#Move boost into install directory
mv boost_1_55_0/* src
rm -rf boost_1_55_0

patch src/boost/mpl/aux_/config/adl.hpp ./configs/1.55/adl.hpp.patch

#Configure and compile
cd src 
./bootstrap.sh \
	--prefix="$ISSM_DIR/externalpackages/boost/install" \
	--with-python=python2.7 \
	--with-python-root="$ISSM_DIR/externalpackages/python/install" 

#Compile boost
./bjam install

#put bjam into install also: 
mkdir ../install/bin
cp bjam ../install/bin
