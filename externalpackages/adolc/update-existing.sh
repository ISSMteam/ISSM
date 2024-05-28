#!/bin/bash
set -eu

#This script is very specific to the larour larsen nightly using adolc. 
#It relies on just an update of adolc, and symlink to an existing adolc 
#repo. 

#Some cleanup
rm -rf install src

#symlink: 
ln -s /proj/ice/larour2/issm-uci/trunk-jpl/externalpackages/adolc/adolc_issm ./src

#update and compile
cd src
git pull origin
git checkout ampi
#git reset --hard b254b2a001a1b7a024a9184cd087ae06eb975cad

autoreconf -f -i 
./configure --prefix=$ISSM_DIR/externalpackages/adolc/install \
	--libdir=$ISSM_DIR/externalpackages/adolc/install/lib 

if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
