#!/bin/bash
set -eu


## Constants
#
VER=2.2

PREFIX="${ISSM_DIR}/externalpackages/chaco/install" # Set to location where external package should be installed

## Environment
#
export CFLAGS="-Wno-error=implicit-function-declaration"

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/Chaco-${VER}.tar.gz" "Chaco-${VER}.tar.gz"
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/docs/chaco_guide.pdf" "chaco_guide.pdf"

# Unpack source
tar -xvzf Chaco-${VER}.tar.gz

# Move source to 'src' directory
mv Chaco-${VER}/* src
rm -rf Chaco-${VER}

# Apply patches
patch -R -p0 < chaco.patch # Written by diff -rc src ~/Libs/Chaco-${VER} > chaco.patch
patch src/code/Makefile < patches/Makefile.patch

# Compile
cd src/code
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make chacominusblas.a

# Clean up objects (but not library or executable)
make clean
cd ../..

# Install
cp -p src/exec/README ${PREFIX}
cp -p src/exec/User_Params ${PREFIX}
cp -p src/exec/*.coords ${PREFIX}
cp -p src/exec/*.graph ${PREFIX}
mkdir ${PREFIX}/include
cp -p src/code/main/defs.h ${PREFIX}/include/defs.h
cp -p src/code/main/params.h ${PREFIX}/include/params.h
cp -p chaco.h ${PREFIX}/include/chaco.h
mkdir ${PREFIX}/lib
mv src/code/chaco.a ${PREFIX}/lib/libchaco.a
mv src/code/chacominusblas.a ${PREFIX}/lib/libchacominusblas.a
mkdir ${PREFIX}/exec
mv src/exec/chaco ${PREFIX}/exec
