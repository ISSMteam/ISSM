#!/bin/bash
set -eu

#Some cleanup
rm -rf install src

#Almost the same as adjoinable mpi but some adaptations to CoDiPack 
svn co https://github.com/michel2323/AdjointMPI.git/branches/doubleFixMaster src
mkdir install

cd src
autoreconf -fi
 ./configure \
	 --prefix="$ISSM_DIR/externalpackages/adjointmpi/install" \
	 --with-mpi-root="/nasa/sgi/mpt/2.15r20/" \
	 CXXFLAGS="-O2 -fPIC" CFLAGS="-O2 -fPIC" 

#Compile adjoinablempi 
make clean
if [ $# -eq 0 ]; then
	make 
else
	make -j $1
fi
make install
