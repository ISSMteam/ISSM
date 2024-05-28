#!/bin/bash
#set -eu # Do not `run set -eu` because it causes some targets to fail


## Constants
#
VER="1_73_0"

PREFIX="${ISSM_DIR}/externalpackages/boost/install" # Set to location where external package should be installed

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/boost_${VER}.tar.gz" "boost_${VER}.tar.gz"

# Unpack source
tar -zxvf boost_${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Move source into 'src' directory
mv boost_${VER}/* src
rm -rf boost_${VER}

# Configure
cd src
./bootstrap.sh \
	--prefix=${PREFIX} \
	--with-python=python2.7

# Modify project config to enable MPI
printf "\n# Enable MPI\nusing mpi ;\n" >> project-config.jam

# Compile and install
./b2 install link=static runtime-link=static

# Remove any dynamic libraries that may have been compiled
#
# TODO: Reconfigure so that dynamic libraries are not compiled at all
#
rm -f $(find ${PREFIX}/lib -name *.so*)
