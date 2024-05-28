#!/bin/bash
set -eu

# NOTE: Stop after bootstrap phase, and run,
#
#	bjam --debug-configuration
#
# to figure out which paths Boost is using to include Python. Make sure 
# every one of these paths is covered by Python. If not, make symlinks in 
# externalpackages/python to what Boost is expecting (assumes that you are 
# using copy of Python installed as an external package of ISSM). There is 
# NO WAY to get the Boost library to include Python support under this 
# configuration without doing this
#

export BOOST_ROOT="${ISSM_DIR}/externalpackages/boost"
export CXXFLAGS="-D_INTEL_LINUX_ -std=c++11"
export CFLAGS="-D_INTEL_LINUX_"

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
# Copy customized source and configuration files to 'src' directory
cp configs/1.55/boost/multi_index/ordered_index.hpp src/boost/multi_index

#Configure and compile
cd src
./bootstrap.sh \
	--prefix=${BOOST_ROOT}/install \
	--with-python=python2.7

#Compile boost
./bjam install

#put bjam into install also: 
mkdir ../install/bin
cp bjam ../install/bin
