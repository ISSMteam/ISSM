#!/bin/bash
#set -eu # Do not `run set -eu` because it causes some targets to fail


## Constants
#
VER="1_55_0"

## Environment
#
export BOOST_ROOT="${ISSM_DIR}/externalpackages/boost"
export CXXFLAGS='-std=c++98' # Setting CXXFLAGS to deal with C++11 incompatibility with MATLAB's Boost

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/boost_${VER}.tar.gz" "boost_${VER}.tar.gz"

# Unpack source
tar -zxvf boost_${VER}.tar.gz

# Cleanup
rm -rf install src
mkdir install src

# Move source into 'src' directory
mv boost_${VER}/* src/
rm -rf boost_${VER}

# Copy customized source and configuration files to 'src' directory
cp configs/1.55/boost/multi_index/ordered_index.hpp src/boost/multi_index

# Configure
cd src
./bootstrap.sh \
	--prefix=${BOOST_ROOT}/install \
	--with-python=python2.7

# Modify project config to enable MPI
printf "\n# Enable MPI\nusing mpi ;\n" >> project-config.jam

# Compile and install
./bjam link=static install

# Copy binary to install directory
mkdir ${BOOST_ROOT}/install/bin
cp bjam ${BOOST_ROOT}/install/bin
