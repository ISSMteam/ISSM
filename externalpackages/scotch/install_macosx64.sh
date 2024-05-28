#!/bin/bash
set -eu

# Some cleanup
rm -rf scotch_5.1
rm -rf src 
rm -rf install 

# Create src and install directories
mkdir src install 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/scotch_5.1.6.tar.gz' 'scotch_5.1.6.tar.gz'
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/docs/ptscotch_user5.1.pdf' 'ptscotch_user5.1.pdf'
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/docs/scotch_user5.1.pdf' 'scotch_user5.1.pdf'

# Untar 
tar -xvzf scotch_5.1.6.tar.gz

# Apply patches (all at once, since many)
# (written by diff -rc old_src new_src > scotch.patch)
patch -p0 < scotch.patch

# Move scotch to src directory
mv scotch_5.1/* src
rm -rf scotch_5.1

# Build scotch
cp -p Makefile.inc.mac src/src/Makefile.inc
cp -p gmap_mex.c src/src/scotch
cp -p gmapx.c src/src/scotch
cd src/src
# For stand-alone scotch modules:
make scotch
make clean
# For no-file-io scotch modules:
make nfioscotch
# Clean up
make clean
cd ../..

# Populate install directory
cp -pr src/grf install
cp -pr src/tgt install
cp -pr src/doc install
cp -pr src/man install
mkdir install/include
cp -p src/src/libscotch/module.h install/include/scotch_module.h
cp -p src/src/libscotch/common.h install/include/scotch_common.h
cp -p src/include/scotch.h install/include/
cp -p src/src/scotch/gmap.h install/include/scotch_gmap.h
mkdir install/lib
mv src/lib/* install/lib
mkdir install/bin
mv src/bin/* install/bin
#cp -p gmap.m install/bin
#cp -p gpart.m install/bin
