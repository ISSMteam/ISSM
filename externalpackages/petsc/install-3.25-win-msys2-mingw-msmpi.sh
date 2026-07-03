#!/bin/bash
set -u # NOTE: Do not set -e as it will cause this script to fail when there are errors in underlying Python scripts

# NOTE:
# - You must install various needed packages with,
#
#		pacman -S mingw-w64-x86_64-msmpi mingw-w64-x86_64-toolchain python
#
# - You must use MSYS2 MinGW 64-bit version of cmake to be able to install 
#	external packages correctly,
#
#		pacman -R mingw-w64-x86_64-cmake
#
# Sources:
# - https://gitlab.com/petsc/petsc/-/issues/820#note_487483240
#

## Constants
#
VER="3.25.3"

PETSC_ARCH="arch-mswin-c-opt"
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

# Patch source
sed -i 's|(MAKEFLAGS)|(MAKEFLAGS:w=)|' ${PETSC_DIR}/makefile ${PETSC_DIR}/lib/petsc/conf/rules # Fix for issue with GNUMake 4.4.1 (https://gitlab.com/petsc/petsc/-/merge_requests/6140)

# Configure
#
cd ${PETSC_DIR}
/usr/bin/python ./configure \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--PETSC_ARCH="${PETSC_ARCH}" \
	--with-shared-libraries=0 \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-mpiexec="${MPIEXEC_DIR}/mpiexec.exe" \
	--download-fblaslapack=1 \
	--download-metis=1 \
	--download-parmetis=1 \
	--download-scalapack=1 \
	--download-mumps=1

# Compile and install
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH="${PETSC_ARCH}" all
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH="${PETSC_ARCH}" install
