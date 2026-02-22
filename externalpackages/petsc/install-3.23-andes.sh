#!/bin/bash
set -eu

## Constants
VER="3.23.6"

PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src" # DO NOT CHANGE THIS
PREFIX="${ISSM_DIR}/externalpackages/petsc/install" # Set to location where external package should be installed

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${VER}.tar.gz" "petsc-${VER}.tar.gz"

# Unpack source
tar -zxvf petsc-${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} ${PETSC_DIR}
mkdir -p ${PETSC_DIR}

# Move source to $PETSC_DIR
mv petsc-${VER}/* ${PETSC_DIR}
rm -rf petsc-${VER}

#CC=icx CXX=icx FC=ifx \

# Configure
cd ${PETSC_DIR}
./configure \
	CC=icx CXX=icpx FC=ifx \
	COPTFLAGS="-g -O3" CXXOPTFLAGS="-g -O3" FOPTFLAGS="-g -O3" \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--with-make-np=20 \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-pic=1 \
	--with-blas-lapack-dir=$MKL_ROOT \
	--download-metis=1 \
	--download-mpich=1 \
	--download-parmetis=1 \
	--download-scalapack=1 \
	--download-mumps=1

# Compile and install
make
make install
