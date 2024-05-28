#!/bin/bash
set -eu
#you need hdf5 compiled

#Some cleanup
rm -rf install 
mkdir install 

git clone https://github.com/nco/nco.git install
cd install
git checkout 4.7.9

#Configure and compile
./configure \
 --disable-doc \
 --prefix="$ISSM_DIR/externalpackages/nco/install" 

if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
