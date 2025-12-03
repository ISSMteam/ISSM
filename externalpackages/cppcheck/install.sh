#!/bin/bash
set -eu


## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/cppcheck/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX} src

# Download source
git clone git@github.com:danmar/cppcheck.git src

# Compile and install
cd src
if [ $# -eq 0 ]; then
	make MATCHCOMPILER=yes FILESDIR="${PREFIX}" HAVE_RULES=yes CXXFLAGS="-O2 -DNDEBUG -Wall -Wno-sign-compare -Wno-unused-function"
else
	make -j $1 MATCHCOMPILER=yes FILESDIR="${PREFIX}" HAVE_RULES=yes CXXFLAGS="-O2 -DNDEBUG -Wall -Wno-sign-compare -Wno-unused-function"
fi
