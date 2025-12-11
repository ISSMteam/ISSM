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
#
# NOTE:
# - Cannot use --with-fpic option when compiling static libs,
#
#		Cannot determine compiler PIC flags if shared libraries is turned off
#		Either run using --with-shared-libraries or --with-pic=0 and supply the
#		compiler PIC flag via CFLAGS, CXXXFLAGS, and FCFLAGS
#
cd ${PETSC_DIR}
./config/configure.py \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--with-shared-libraries=0 \
	--CFLAGS="-fPIC" \
	--CXXFLAGS="-fPIC" \
	--FFLAGS="-fPIC" \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--download-fblaslapack=1 \
	--download-mpich=1 \
	--download-metis=1 \
	--download-parmetis=1 \
	--download-scalapack=1 \
	--download-mumps=1 \
	--download-zlib=1 \
	--download-hdf5=1

# Compile and install
make
make install
