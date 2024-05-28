#!/bin/bash
set -eu

##
# Source code and installation guidelines available at,
#
#	https://github.com/nakatamaho/mplapack
#
# Note that original website (for MPACK) at,
#
#	http://mplapack.sourceforge.net
#
# is not maintained.
#

## Constants
#
VER="0.9.3"

PREFIX="${ISSM_DIR}/externalpackages/mplapack/install" # Set to location where external package should be installed

## Environment
#
export CC=gcc
export CXX=g++
export FC=gfortran
export F77=gfortran

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/mplapack-${VER}.tar.gz" "mplapack-${VER}.tar.gz"

# Unpack source
tar -zxvf mplapack-${VER}.tar.gz

# Move source to 'src' directory
mv mplapack-${VER}/* src
rm -rf mplapack-${VER}

# Configure
cd src
pushd mplapack/debug
bash gen.Makefile.am.sh
popd
autoreconf -ivf
./configure \
	--prefix="${PREFIX}" \
	--disable-static \
	--enable-fast-install=yes \
	--disable-dependency-tracking \
	--enable-debug=yes \
	--enable-fortranwrapper=yes \
	--enable-gmp=yes \
	--enable-mpfr=yes \
	--enable-qd=yes \
	--enable-dd=yes \
	--enable-double=yes \
	--enable-_Float128=yes \
	--enable-_Float64x=yes \
	--with-external-blas="-L${BLAS_ROOT}/lib -lfblas" \
	--with-external-lapack="-L${LAPACK_ROOT}/lib -lflapack"

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
