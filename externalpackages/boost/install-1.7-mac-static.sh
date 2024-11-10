#!/bin/bash
#set -eu # Do not `run set -eu` because it causes some targets to fail


## Constants
#
VER="1.73.0"

PREFIX="${ISSM_DIR}/externalpackages/boost/install" # Set to location where external package should be installed

## Environment
#
export LDFLAGS="-Wl,-headerpad_max_install_names"

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
rm -f $(find ${PREFIX}/lib -name *.dylib)
