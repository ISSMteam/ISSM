#!/bin/bash
set -eu

#Note of caution:  stop after boostrap phase, and run 
#bjam --debug-configuration, to figure out which paths boost is using to include 
#python. make sure everyone of these paths is covered by python. If not, just make 
#symlinks in externalpackages/python to what boost is expecting. Ther is NO WAY 
#to get the boost library to include python support without doing that. 

## Constants
#
VER="1.55.0"

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

#Setting CXXFLAGS to deal with C++11 incompatibility with MATLAB's Boost
#export PATH="/usr/bin":$PATH
export CXXFLAGS='-std=c++98'
export CC=mpicc
export CXX=mpicxx

#Configure and compile
cd src 
./bootstrap.sh \
	--prefix="$ISSM_DIR/externalpackages/boost/install" \
	--with-python=python3.2 

#Compile boost
#./bjam install
./bjam install

#put bjam into install also: 
mkdir ../install/bin
cp bjam ../install/bin

