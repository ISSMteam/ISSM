#!/bin/bash
set -eu


# WARNING: Make sure you have the right MPI

## Constants
#
VER="3.14.0"

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

#configure
cd ${PETSC_DIR}
./config/configure.py \
	COPTFLAGS="-g -O3" CXXOPTFLAGS="-g -O3" FOPTFLAGS="-g -O3" \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--with-blas-lapack-dir="/sopt/INTEL/compilers_and_libraries_2018.3.222/linux/mkl/" \
	--with-mpi-dir="/sopt/OpenMPI/3.1.2/intel-2018.3-slim/" \
	--known-mpi-shared-libraries=1 \
	--known-64-bit-blas-indices \
	--known-mpi-long-double=1 \
	--known-mpi-int64_t=1 \
	--known-mpi-c-double-complex=1 \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-batch=1  \
	--with-shared-libraries=1 \
	--download-metis=1 \
	--download-parmetis=1 \
	--download-scalapack=1 \
	--download-mumps=1 

# Compile and install
make
make install
