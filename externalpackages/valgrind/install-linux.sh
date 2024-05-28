#!/bin/bash
set -eu


## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/valgrind/install" # Set to location where external package should be installed

# Clean up
rm -rf ${PREFIX} src

# Download development version (the current release never supports the latest 
# OS X releases)
git clone git://sourceware.org/git/valgrind.git src

# Configure
cd src
./autogen.sh
./configure \
	--prefix="${PREFIX}" \
	--enable-only64bit

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make install -j $1
fi
