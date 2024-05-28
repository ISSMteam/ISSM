#!/bin/bash
set -eu

# Cleanup from previous installation
rm -rf install 

#Download development version
git clone https://bitbucket.org/mituq/muq2 install

#configure
cd install
cmake -DCMAKE_INSTALL_PREFIX="${ISSM_DIR}/externalpackages/muq/install"

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
