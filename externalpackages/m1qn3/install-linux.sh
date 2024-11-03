#!/bin/bash
set -eu

## Constants
#
VER=3.3

PREFIX="${ISSM_DIR}/externalpackages/m1qn3/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
#$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/m1qn3-${VER}-distrib.tgz" "m1qn3-${VER}-distrib.tgz"

# Unpack source
tar -xzf m1qn3-${VER}-distrib.tgz

# Move source to 'src' directory
mv m1qn3-${VER}-distrib/* src
rm -rf m1qn3-${VER}-distrib

# Apply patches to src
patch src/src/m1qn3.f patch/m1qn3.f.patch
patch src/blas/ddot.f patch/ddot.f.patch

# Compile and install
if which mpifort 2>/dev/null; then
	export FC="mpifort"
else
	if which ifort 2>/dev/null; then
		export FC="ifort"
	else
		export FC="gfortran"
	fi
fi
export FFLAGS="-fPIC" # -O2 is default

cd src/src
cp ../../configs/makefile .
cp ../../configs/linux/configure.make .
export LIBNAME="m1qn3"
make shared
cp lib${LIBNAME}.* ${PREFIX}

cd ../blas
cp ../../configs/makefile .
cp ../../configs/linux/configure.make .
export LIBNAME="ddot"
make shared
cp lib${LIBNAME}.* ${PREFIX}
