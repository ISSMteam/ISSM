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

# Modify source so that Python 3 can be used to compile PETSc
sed -i'' 's|#!/usr/bin/env python|#!/usr/bin/env python3|g' ${PETSC_DIR}/config/configure.py

# Modify source so that Python >= 3.9 can be used to compile PETSc
sed -i'' 's|thread.isAlive|thread.is_alive|g' ${PETSC_DIR}/config/BuildSystem/script.py

# Configure
#
# NOTE:
# - Added -Wno-error=implicit-function-declaration to CFLAGS for Clang >= 12. 
#	(may need to remove it for earlier versions not using the C99 standard).
#
cd ${PETSC_DIR}
./config/configure.py \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--CFLAGS="-Wno-error=implicit-function-declaration" \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-pic=1 \
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
