#!/bin/bash
set -eu

#Erase install
rm -rf install
mkdir install

#Download from GitHub server
#GIT must be installed first. See $ISSM_DIR/externalpackages/git
git clone https://github.com/labmec/neopz.git

#use these 3 lines if one would like to change to a stable version
cd neopz
git checkout 5b1d4fa3cf61dcc500742b8cfdfb01d86ec724b3
cd ..

#Untar and set src directory
mv neopz/ install/

#Set neopz CMake variables
#CMake must be installed first. See $ISSM_DIR/externalpackages/git
export PROJECT_SOURCE_DIR=$ISSM_DIR/externalpackages/neopz/install/neopz
export PROJECT_BINARY_DIR=$ISSM_DIR/externalpackages/neopz/install/

#Configure neopz using cmake
cd $PROJECT_SOURCE_DIR
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PROJECT_BINARY_DIR -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-g -O3"

cd $PROJECT_SOURCE_DIR
make
make install

cd $PROJECT_BINARY_DIR/pzlib
mv lib ../
mv include ../
cd ..
rm -rf pzlib
cd ..
