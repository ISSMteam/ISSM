#!/bin/bash
set -eu


## Constants
#
VER="3.14.6"

PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src" # DO NOT CHANGE THIS
PREFIX="${ISSM_DIR}/externalpackages/petsc/install" # Set to location where external package should be installed

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${VER}.tar.gz" "petsc-${VER}.tar.gz"

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
./config/configure.py \
	COPTFLAGS="-g -O2" CXXOPTFLAGS="-g -O2" FOPTFLAGS="-g -O2" \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-pic=1 \
	--with-blas-lapack-dir="/software/sse/easybuild/prefix/software/imkl/2018.1.163-iimpi-2018a/mkl" \
	--with-cc="/software/sse/easybuild/prefix/software/impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28/bin64/mpicc" \
	--with-cxx="/software/sse/easybuild/prefix/software/impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28/bin64/mpicxx" \
	--with-fc="/software/sse/easybuild/prefix/software/impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28/bin64/mpiifort" \
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
