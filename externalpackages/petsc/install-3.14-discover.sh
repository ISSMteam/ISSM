#!/bin/bash
set -eu


## Constants
#
VER="3.14.6"

PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src" # DO NOT CHANGE THIS
PREFIX="${ISSM_DIR}/externalpackages/petsc/install" # Set to location where external package should be installed

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-lite-${VER}.tar.gz" "petsc-${VER}.tar.gz"

# Unpack source
tar -zxvf petsc-${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} ${PETSC_DIR}
mkdir -p ${PETSC_DIR}

# Move source to $PETSC_DIR
mv petsc-${VER}/* ${PETSC_DIR}
rm -rf petsc-${VER}

# Configure
cd ${PETSC_DIR}
./configure \
	COPTFLAGS="-g -O3" CXXOPTFLAGS="-g -O3" FOPTFLAGS="-g -O3" \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-pic=1 \
	--with-blas-lapack-dir="/usr/local/intel/oneapi/2021/mkl/2021.4.0/" \
	--with-cc="/usr/local/intel/oneapi/2021/mpi/2021.4.0/bin/mpicc" \
	--with-cxx="/usr/local/intel/oneapi/2021/mpi/2021.4.0/bin/mpicxx" \
	--with-fc="/usr/local/intel/oneapi/2021/mpi/2021.4.0/bin/mpif90" \
	--known-mpi-shared-libraries=1 \
	--known-64-bit-blas-indices \
	--known-mpi-long-double=1 \
	--known-mpi-int64_t=1 \
	--known-mpi-c-double-complex=1 \
	--with-shared-libraries=1 \
	--download-metis=1 \
	--download-parmetis=1 \
	--download-scalapack=1 \
	--download-mumps=1

# Compile and install
make
make install
