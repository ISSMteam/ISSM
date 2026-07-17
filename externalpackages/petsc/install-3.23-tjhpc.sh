#!/bin/bash
set -eu


## Constants
#
VER="3.23.6"

PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src" # DO NOT CHANGE THIS
PREFIX="${ISSM_DIR}/externalpackages/petsc/install" # Set to location where external package should be installed

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${VER}.tar.gz" "petsc-${VER}.tar.gz"

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
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--with-make-np=20 \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-pic=1 \
	--download-fblaslapack="${ISSM_DIR}/externalpackages/petsc/downloads/petsc-pkg-fblaslapack-e8a03f57d64c.tar.gz" \
	--download-metis="${ISSM_DIR}/externalpackages/petsc/downloads/petsc-pkg-metis-69fb26dd0428.tar.gz" \
	--download-mpich="${ISSM_DIR}/externalpackages/petsc/downloads/mpich-4.3.0.tar.gz" \
	--download-mumps="${ISSM_DIR}/externalpackages/petsc/downloads/MUMPS_5.7.3.tar.gz" \
	--download-parmetis="${ISSM_DIR}/externalpackages/petsc/downloads/petsc-pkg-parmetis-f5e3aab04fd5.tar.gz" \
	--download-scalapack="${ISSM_DIR}/externalpackages/petsc/downloads/scalapack-0e8767285b7a201c7b1ff34d2c2bb009534145df.tar.gz" \
	--download-zlib="${ISSM_DIR}/externalpackages/petsc/downloads/zlib-1.3.1.tar.gz"

# Compile and install
make
make install
