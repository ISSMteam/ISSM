#!/bin/bash
#set -eu # Do not `run set -eu` because it causes some targets to fail


## Constants
#
VER="1_73_0"

PREFIX="${ISSM_DIR}/externalpackages/boost/install" # Set to location where external package should be installed

## Environment
#
export LDFLAGS="-Wl,-headerpad_max_install_names"

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
./b2 install link=shared runtime-link=shared

# Set install_name for all shared libraries
#
# NOTE: The install_name_tool prints an error message on some installations, 
#		but this is not a shell error, and it does not seem to affect the 
#		Boost libraries called by ISSM. For now, we are simply redirecting the 
#		error to null.
#
# TODO:
# - Modify the source to apply absolute paths to the library ids so that 
#	patching it after the fact with install_name_tool is not necessary.
#
cd ${PREFIX}/lib
for name in *.dylib; do
	install_name_tool -id ${PREFIX}/lib/${name} ${name} 2>/dev/null
done

if [ "${VER}" == "1_79_0" ]; then
	## Patch install names for certain libraries
	#
	# TODO: Figure out how to reconfigure source to apply these install names at compile time
	#
	install_name_tool -change @rpath/libboost_atomic.dylib ${PREFIX}/lib/libboost_atomic.dylib libboost_filesystem.dylib
fi
