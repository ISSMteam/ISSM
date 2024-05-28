#!/bin/bash

cd adolc_issm
git pull

autoreconf -f -i 
./configure --prefix=$ISSM_DIR/externalpackages/adolc/install 

if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
