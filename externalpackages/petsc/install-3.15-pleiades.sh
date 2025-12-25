#!/bin/bash
set -eu


## Constants
#
VER="3.15.0"

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
#
# NOTE:
# - Based on /nasa/petsc/3.7.5/intel_mpt/lib/petsc/conf/petscvariables (look 
#	for CONFIGURE_OPTIONS)
#
cd ${PETSC_DIR}
./configure \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--with-blas-lapack-dir="/nasa/intel/Compiler/2018.3.222/compilers_and_libraries_2018.3.222/linux/mkl" \
	--with-scalapack-include=/nasa/intel/Compiler/2018.3.222/mkl/include \
	--with-scalapack-lib="/nasa/intel/Compiler/2018.3.222/mkl/lib/intel64/libmkl_scalapack_lp64.so /nasa/intel/Compiler/2018.3.222/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so" \
	--with-make-np=10 \
	--known-mpi-shared-libraries=1 \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-batch=1  \
	--with-shared-libraries=1 \
	--download-metis=1 \
	--download-parmetis=1 \
	--download-mumps=1 \
	--download-scalapack=0

# Compile and install
make
make install
