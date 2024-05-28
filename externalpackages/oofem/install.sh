#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf oofem-2.0
mkdir src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/oofem-2.0.tar.gz' 'oofem-2.0.tar.gz'

#Untar 
tar -zxvf  oofem-2.0.tar.gz

#Move oofem into src directory
mv oofem-2.0/* src
rm -rf oofem-2.0

# currently a basic serial configuration, see http://www.oofem.org/wiki/doku.php?id=installation 
# for details on enabling the IML++ and PETSc libraries, and for configuring the parallel version
cd src
./configure OOFEM_TARGET=oofem-2.0 --enable-dss

#Compile oofem 
cd targets/oofem-2.0
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi

# build verification tools
cd ../../tools && make all

# testing solver
cd ../targets/oofem-2.0 && make tests && less ./test_results
