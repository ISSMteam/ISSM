#!/bin/bash
set -eu

#Erase install
rm -rf install  src 

#Download
git clone git://github.com/danmar/cppcheck.git src

#compile
cd src
if [ $# -eq 0 ]; then
	make 
else 
	make -j $1
fi  
make install PREFIX=$ISSM_DIR/externalpackages/cppcheck/install
cd ..
