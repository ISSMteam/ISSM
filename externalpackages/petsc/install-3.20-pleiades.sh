#!/bin/bash
set -eu


## Constants
#
VER="3.20.1"

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

# Configure
#
# NOTE:
# - Options from,
#
# 		cat /nasa/petsc/3.14.5_toss3/lib/petsc/conf/petscvariables | grep CONF
#
cd ${PETSC_DIR}
./configure \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--CFLAGS="-g -O3" \
	--CXXFLAGS="-g -O3" \
	--FFLAGS="-g -O3" \
	--with-make-np=10 \
	--with-batch=1 \
	--with-pic=1 \
	--with-shared-libraries=1 \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--download-fblaslapack=1 \
	--download-scalapack=1 \
	--download-metis=1 \
	--download-parmetis=1 \
	--download-mumps=1 

# Compile and install
make
make install
