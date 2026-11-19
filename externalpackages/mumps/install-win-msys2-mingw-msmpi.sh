#!/bin/bash
set -eu

## Constants
#
VER="5.9.1"

PREFIX="${ISSM_DIR}/externalpackages/mumps/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} ${PREFIX}/include ${PREFIX}/lib src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "http://mumps-solver.org/MUMPS_${VER}.tar.gz" "MUMPS_${VER}.tar.gz"

# Unpack source
tar -zxvf MUMPS_${VER}.tar.gz

# Move source to 'src' directory
mv MUMPS_${VER}/* src
rm -rf MUMPS_${VER}

# Copy Makefile to 'src' directory
cp configs/${VER}/Make.inc/Makefile.WIN-MSYS2-MINGW.PAR src/Makefile.inc

# Patch source where needed
patch src/src/cana_aux_par.F < configs/${VER}/src/cana_aux_par.F.patch
patch src/src/dana_aux_par.F < configs/${VER}/src/dana_aux_par.F.patch
patch src/src/sana_aux_par.F < configs/${VER}/src/sana_aux_par.F.patch
patch src/src/zana_aux_par.F < configs/${VER}/src/zana_aux_par.F.patch

# Compile and install
cd src
make all
cd ..
cp src/include/* ${PREFIX}/include
cp src/lib/* ${PREFIX}/lib
