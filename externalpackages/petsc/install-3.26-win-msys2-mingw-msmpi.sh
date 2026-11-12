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
VER="3.26.1"

PETSC_ARCH_INIT="arch-mswin-init"
PETSC_ARCH_RECONF="arch-mswin-reconf"
PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src" # DO NOT CHANGE THIS
PREFIX="${ISSM_DIR}/externalpackages/petsc/install" # Set to location where external package should be installed (UNIX)
PREFIX_WIN="${ISSM_DIR_WIN}/externalpackages/petsc/install" # Set to location where external package should be installed (POSIX)

## Environment
#
export TMPDIR=/tmp # Fixes permission issue when running CMake when user is not an administrator
export TMP=/tmp
export TEMP=/tmp

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${VER}.tar.gz" "petsc-${VER}.tar.gz"

# Unpack source
tar -zxvf petsc-${VER}.tar.gz

# Cleanup
rm -rf ${PREFIX} ${PETSC_DIR}
mkdir -p ${PREFIX} ${PETSC_DIR}

# Move source to $PETSC_DIR
mv petsc-${VER}/* ${PETSC_DIR}
rm -rf petsc-${VER}

# Patch source
sed -i 's|(MAKEFLAGS)|(MAKEFLAGS:w=)|' ${PETSC_DIR}/makefile ${PETSC_DIR}/lib/petsc/conf/rules # Fix for issue with GNUMake 4.4.1 (https://gitlab.com/petsc/petsc/-/merge_requests/6140)

# Configure
#
# NOTE: --prefix="${PREFIX_WIN}" is necessary to, at least, configure and 
#       compile ScaLAPACK under UCRT.
#
cd ${PETSC_DIR}
echo "Configuring with Windows-style prefix before compiling..."
/usr/bin/python ./configure \
    --prefix="${PREFIX_WIN}" \
    --PETSC_ARCH="${PETSC_ARCH_INIT}" \
    --with-shared-libraries=0 \
    --with-debugging=0 \
    --with-valgrind=0 \
    --with-x=0 \
    --with-ssl=0 \
    --with-cc=/ucrt64/bin/mpicc.exe \
    --with-cxx=/ucrt64/bin/mpicxx.exe \
    --with-fc=/ucrt64/bin/mpifort.exe \
    --with-mpiexec="${MPIEXEC_DIR}/mpiexec.exe" \
    --download-fblaslapack=1 \
    --download-metis=1 \
    --download-parmetis=1 \
    --download-scalapack=1

# Now reconfigure with the Unix-style path, otherwise the installation of 
# libpetsc fails).
#
echo "Reconfiguring with Unix-style prefix before installation..."
/usr/bin/python ./configure \
    --prefix="${PREFIX}" \
    --PETSC_ARCH="${PETSC_ARCH_RECONF}" \
    --with-shared-libraries=0 \
    --with-debugging=0 \
    --with-valgrind=0 \
    --with-x=0 \
    --with-ssl=0 \
    --with-cc=/ucrt64/bin/mpicc.exe \
    --with-cxx=/ucrt64/bin/mpicxx.exe \
    --with-fc=/ucrt64/bin/mpifort.exe \
    --with-mpiexec="${MPIEXEC_DIR}/mpiexec.exe" \
    --with-blas-lib="-L${PETSC_DIR}/${PETSC_ARCH_INIT}/externalpackages/git.fblaslapack -lfblas" \
    --with-lapack-lib="-L${PETSC_DIR}/${PETSC_ARCH_INIT}/externalpackages/git.fblaslapack -lflapack"

# Compile and install
make PETSC_ARCH="${PETSC_ARCH_RECONF}" all
make PETSC_ARCH="${PETSC_ARCH_RECONF}" install

# Cleanup
rm -rf ${TMPDIR}/petsc
