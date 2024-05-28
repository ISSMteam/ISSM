#!/bin/bash
set -eu

#Some cleanup
rm -rf install valgrind-3.10.0
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/valgrind-3.10.0.tar.bz2' 'valgrind-3.10.0.tar.bz2'

#Untar 
tar -jxvf  valgrind-3.10.0.tar.bz2

#Move valgrind into install directory
mv valgrind-3.10.0/* install
rm -rf valgrind-3.10.0

#configure
cd install
./configure --prefix="$ISSM_DIR/externalpackages/valgrind/install"

#Compile valgrind
make  -j 4
make install

#final thing: if mpi is compiled in, soft link its target to a simpler name
cd lib
ln -s valgrind/libmpi*  ./libmpidebug.so
