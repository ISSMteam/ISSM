#!/bin/bash
set -eu

# Some cleanup
rm -rf shapelib-1.2.10
rm -rf src 
rm -rf install 
mkdir src install 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/shapelib-1.2.10.tar.gz' 'shapelib-1.2.10.tar.gz'
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/shape_eg_data.zip'  'shape_eg_data.zip'
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/docs/shapefile.pdf' 'shapefile.pdf'

# Untar 
tar -xvzf shapelib-1.2.10.tar.gz
unzip shape_eg_data.zip -d shapelib-1.2.10/eg_data

# Move shapelib to src directory
mv shapelib-1.2.10/* src
rm -rf shapelib-1.2.10

# Apply patches (all at once)
# (written by diff -rc old_src new_src > shapelib.patch)
patch -p0 < shapelib.patch

# Build shapelib and run self-tests
cd src
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make test

# Clean up objects (but not library or executable)
#make clean
cd ..

# Populate install directory
cp -p src/README install
cp -p src/*.html install
mkdir install/include
cp -p src/shapefil.h install/include
mkdir install/lib
mv src/libshape.a install/lib/libshape.a
mkdir install/exec
mv src/shpcreate install/exec
mv src/shpadd install/exec
mv src/shpdump install/exec
mv src/shprewind install/exec
mv src/dbfcreate install/exec
mv src/dbfadd install/exec
mv src/dbfdump install/exec
mv src/shputils install/exec
mv src/shptest install/exec
