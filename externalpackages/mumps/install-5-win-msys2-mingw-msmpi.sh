#!/bin/bash
set -eu


# Sources:
# - https://www.scivision.dev/windows-mpi-msys2/
#
# NOTE:
# - Source is pulled from https://github.com/scivision/mumps/releases, which 
#	provides patches to the releases from developer (http://mumps-solver.org/)
#
# TODO:
# - Create install alias in Makefiles
# - Alternatively, use cmake rather than make
#

## Constants
#
VER=5.3.5.2

PREFIX="${ISSM_DIR}/externalpackages/mumps/install"

# Cleanup
rm -rf ${PREFIX} src
mkdir ${PREFIX} src src/lib

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/mumps-${VER}.tar.gz" "mumps-${VER}.tar.gz"

# Unpack source
tar -zxvf mumps-${VER}.tar.gz

# Move source into 'src' directory
mv mumps-${VER}/* src
rm -rf mumps-${VER}

# Copy customized source and configuration files to 'src' directory
cp configs/5.3/win/msys2/mingw64/msmpi/Makefile.inc src
cp configs/5.3/win/msys2/mingw64/msmpi/libseq/Makefile src/libseq
cp configs/5.3/win/msys2/mingw64/msmpi/PORD/lib/Makefile src/PORD/lib
cp configs/5.3/win/msys2/mingw64/msmpi/src/Makefile src/src

# Compile
cd src
if [ $# -eq 0 ]; then
	make all
else
	make -j $1 all
fi

# Install
mkdir ${PREFIX}/lib
cp lib/lib* ${PREFIX}/lib
cp libseq/lib* ${PREFIX}/lib
mkdir ${PREFIX}/include
cp include/* ${PREFIX}/include

# Create link to lib directory (PETSc, by default, looks for libraries in 
# lib64/ if it detects that 64-bit integers are being used).
cd ${PREFIX}
ln -s lib lib64
