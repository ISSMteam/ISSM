#!/bin/bash
set -eu

#Some cleanup
rm -rf Dakota
rm -rf src 
rm -rf install 
mkdir src install 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/Dakota_5_2.src.tar.gz' 'Dakota_5_2.src.tar.gz'

#Untar 
tar -zxvf Dakota_5_2.src.tar.gz

#Move Dakota to src directory
mv Dakota/* src
rm -rf Dakota

#Apply patches
patch src/src/ParallelLibrary.C configs/5.2/ParallelLibrary.C.patch
patch src/src/ParallelLibrary.H configs/5.2/ParallelLibrary.H.patch
#patch src/src/NIDRProblemDescDB.C configs/5.2/NIDRProblemDescDB.C.patch    #  causes segfault in 5.2
patch src/src/NonDSampling.C configs/5.2/NonDSampling.C.patch
patch src/src/NonDLocalReliability.C configs/5.2/NonDLocalReliability.C.patch
patch src/src/NonDUnilevelRBDO.C configs/5.2/NonDUnilevelRBDO.C.patch    #  source not even used?
patch src/packages/pecos/src/pecos_global_defs.hpp configs/5.2/pecos_global_defs.hpp.patch

#Configure dakota
cd src
./configure \
	--prefix="$ISSM_DIR/externalpackages/dakota/install" \
	--without-graphics  \
	--with-pic \
	--disable-mpi \
	--with-plugin \
	--with-blas="-L$ISSM_DIR/externalpackages/petsc/install/lib -lfblas " \
	--with-lapack="-L$ISSM_DIR/externalpackages/petsc/install/lib -lflapack "
cd ..

#Compile and install dakota
cd src 
if [ $# -eq 0 ];
then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
cd ..
