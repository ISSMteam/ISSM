#!/bin/bash
set -u # NOTE: Do not set -e as it will cause this script to fail when there are errors in underlying Python scripts

# NOTE:
# - You must install various needed packages with,
#
#		pacman -S mingw-w64-x86_64-toolchain python
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
VER="3.14.6"

MAKEFILE_GENERATOR='-G "MSYS Makefiles"'
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
sed -i -e 's|#!/usr/bin/env python|#!/usr/bin/python3|' ${PETSC_DIR}/config/configure.py
sed -i 's|  chkusingwindowspython()|#  chkusingwindowspython()|' ${PETSC_DIR}/config/configure.py
sed -i 's|(MAKEFLAGS)|(MAKEFLAGS:w=)|' ${PETSC_DIR}/makefile ${PETSC_DIR}/lib/petsc/conf/rules # Fix for issue with GNUMake 4.4.1 (https://gitlab.com/petsc/petsc/-/merge_requests/6140)

# Configure
# - Cannot use --with-fpic option when compiling static libs,
#
#		Cannot determine compiler PIC flags if shared libraries is turned off
#		Either run using --with-shared-libraries or --with-pic=0 and supply the
#		compiler PIC flag via CFLAGS, CXXXFLAGS, and FCFLAGS
#
# - Added -fallow-argument-mismatch option to FFLAGS in order to clear "Error: 
#	Rank mismatch between actual argument at [...]"
# - Added -fallow-invalid-boz option to FFLAGS in order to clear "Error: BOZ 
#	literal constant at [...]"
# - Argument to --with-mpi-include must be a list or it gets expanded 
#	incorrectly
#
cd ${PETSC_DIR}
./config/configure.py \
	--prefix="${PREFIX}" \
	--PETSC_DIR="${PETSC_DIR}" \
	--PETSC_ARCH="${PETSC_ARCH}" \
	--CFLAGS="-fPIC -Wl,-static -Wno-error=implicit-function-declaration -Wno-error=implicit-int" \
	--CXXFLAGS="-fPIC -Wl,-static" \
	--FFLAGS="-fPIC -Wl,-static -fallow-argument-mismatch -fallow-invalid-boz" \
	--with-shared-libraries=0 \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-proc-filesystem=0 \
	--with-mpiexec="${MPIEXEC_DIR}/mpiexec.exe" \
	--with-mpi-lib="-L${MSMPI_ROOT}/lib -lmsmpi" \
	--with-mpi-include="${MSMPI_ROOT}/include" \
	--download-fblaslapack=1 \
	--download-metis=1 \
	--download-metis-cmake-arguments="${MAKEFILE_GENERATOR}" \
	--download-parmetis=1 \
	--download-parmetis-cmake-arguments="${MAKEFILE_GENERATOR}" \
	--download-scalapack=1 \
	--download-scalapack-cmake-arguments="${MAKEFILE_GENERATOR}" \
	--download-mumps=1

# Compile and install
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH="${PETSC_ARCH}" all
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH="${PETSC_ARCH}" install
