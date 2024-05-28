#!/bin/bash
set -eu

#Erase install
rm -rf source build install

#Download code
hg clone http://mercurial.mcs.anl.gov/ad/ADIC -r 631 source

#Then update to -r 631
cd source 
hg update -r 631

#Configure and Compile
cd source
./autogen.sh
./aclocal
cd ..
mkdir build
cd build
../source/configure \
	--with-rose=$ISSM_DIR/externalpackages/rose/install \
	--with-openanalysis=$ISSM_DIR/externalpackages/openanalysis/openanalysis/x86_64-Linux \
	--with-boost=$ISSM_DIR/externalpackages/boost/install \
	--with-xerces=$ISSM_DIR/externalpackages/xerces/src \
	--with-xaifbooster=$ISSM_DIR/externalpackages/xaifbooster/xaifBooster \
	--with-colpack=$ISSM_DIR/externalpackages/colpack/install\
	--prefix=$ISSM_DIR/externalpackages/adic/install 
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
