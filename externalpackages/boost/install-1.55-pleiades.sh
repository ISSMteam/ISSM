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

## Constants
#
VER="1.55.0"

## Environment
#
export BOOST_ROOT="${ISSM_DIR}/externalpackages/boost"
export CXXFLAGS="-D_INTEL_LINUX_ -std=c++11"
export CFLAGS="-D_INTEL_LINUX_"

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://archives.boost.io/release/${VER}/source/boost_${VER//./_}.tar.gz" "boost_${VER//./_}.tar.gz"

# Unpack source
tar -zxvf boost_${VER//./_}.tar.gz

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Move source into 'src' directory
mv boost_${VER//./_}/* src
rm -rf boost_${VER//./_}

patch src/boost/mpl/aux_/config/adl.hpp ./configs/${VER%.*}/adl.hpp.patch
# Copy customized source and configuration files to 'src' directory
cp configs/${VER%.*}/boost/multi_index/ordered_index.hpp src/boost/multi_index

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
