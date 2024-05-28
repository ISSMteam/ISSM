#!/bin/bash
set -eu

#Some cleanup
rm -rf install src
mkdir install

#Download latest version
git clone https://github.com/doxygen/doxygen.git install

#Configure doxygen
cd install 
mkdir build
cd build
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$ISSM_DIR/externalpackages/doxygen/install" ..
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi

make install
