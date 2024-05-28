#!/bin/bash
set -eu

#Some cleanup
rm -rf Dakota
rm -rf src 
rm -rf install 
mkdir src install 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/Dakota_4_3.src.tar.gz' 'Dakota_4_3.src.tar.gz'

#Untar 
tar -zxvf  Dakota_4_3.src.tar.gz

#Move Dakota to src directory
mv Dakota/* src
rm -rf Dakota

#Apply patches
patch src/src/ParallelLibrary.C configs/4.2/ParallelLibrary.C.patch
patch src/src/ParallelLibrary.H configs/4.2/ParallelLibrary.H.patch
patch src/src/NIDRProblemDescDB.C configs/4.2/NIDRProblemDescDB.C.patch
patch src/src/NonDSampling.C configs/4.2/NonDSampling.C.patch
patch src/src/NonDLocalReliability.C configs/4.2/NonDLocalReliability.C.patch
patch src/src/NonDUnilevelRBDO.C configs/4.2/NonDUnilevelRBDO.C.patch    #  source not even used?
#patch -R src/packages/Pecos/src/LHSDriver.cpp configs/4.2/LHSDriver.cpp.patch

#Configure dakota
cd src
./configure \ 
	--prefix="$ISSM_DIR/externalpackages/dakota/install" \
	--without-graphics  \
	--with-pic \
	--disable-mpi \
	--with-blas=/opt/intel/mkl/9.1.023/lib/64/libmkl.so \
	--with-lapack=/opt/intel/mkl/9.1.023/lib/64/libmkl_lapack.so 
cd ..

#Before compiling, if running on 64 bits, we need to active fPIC compilation. Some packages 
#do not register -fPIC in Dakota, which is a problem. Edit the faulty Makefiles and add the -fPIC 
#flag to the compilation.
cat ./src/methods/NCSUOpt/Makefile | sed 's/FFLAGS = -g -O2/FFLAGS = -g -O2 -fPIC/g' >  temp
mv temp ./src/methods/NCSUOpt/Makefile

cat ./src/methods/acro/packages/pebbl/src/Makefile | sed 's/CXXFLAGS = -O2 -fpermissive/CXXFLAGS = -O2 -fpermissive -fPIC/g' > temp
mv temp ./src/methods/acro/packages/pebbl/src/Makefile

cat ./src/methods/hopspack/src-nappspack/Makefile | sed 's/CXXFLAGS = -g -O2/CXXFLAGS = -g -O2  -fPIC/g' > temp
mv temp ./src/methods/hopspack/src-nappspack/Makefile

cat ./src/methods/hopspack/src-cddlib/Makefile | sed 's/CFLAGS = -g -O2/CFLAGS = -g -O2 -fPIC/g' > temp
mv temp  ./src/methods/hopspack/src-cddlib/Makefile 

cat ./src/methods/hopspack/src-shared/Makefile | sed 's/CFLAGS = -g -O2/CFLAGS = -g -O2 -fPIC/g' > temp
mv temp  ./src/methods/hopspack/src-shared/Makefile 

cat ./src/methods/hopspack/src-shared/Makefile | sed 's/CXXFLAGS = -g -O2/CXXFLAGS = -g -O2  -fPIC/g' > temp
mv temp  ./src/methods/hopspack/src-shared/Makefile 

cat ./src/methods/hopspack/src-conveyor/Makefile | sed 's/CXXFLAGS = -g -O2/CXXFLAGS = -g -O2 -fPIC/g' > temp
mv temp  ./src/methods/hopspack/src-conveyor/Makefile 

cat ./src/methods/hopspack/src-appspack/Makefile | sed 's/CXXFLAGS = -g -O2/CXXFLAGS = -g -O2  -fPIC/g' > temp
mv temp ./src/methods/hopspack/src-appspack/Makefile 

cat ./src/methods/acro/packages/colin/src/Makefile | sed 's/CXXFLAGS = -O2 -fpermissive/CXXFLAGS = -O2 -fpermissive -fPIC/g' > temp
mv temp ./src/methods/acro/packages/colin/src/Makefile

cat ./src/methods/acro/packages/coliny/src/Makefile | sed 's/CXXFLAGS = -O2 -fpermissive/CXXFLAGS = -O2 -fpermissive -fPIC/g' > temp
mv temp ./src/methods/acro/packages/coliny/src/Makefile

cat ./src/methods/acro/packages/tpl/3po/Makefile | sed 's/CFLAGS = -O2/CFLAGS = -O2 -fPIC/g' > temp
mv temp  ./src/methods/acro/packages/tpl/3po/Makefile 

cat ./src/methods/acro/packages/tpl/3po/Makefile | sed 's/CXXFLAGS = -O2 -fpermissive/CFLAGS = -O2 -fpermissive -fPIC/g' > temp
mv temp  ./src/methods/acro/packages/tpl/3po/Makefile 

cat ./src/packages/ampl/Makefile | sed 's/CFLAGS = -g -O2/CFLAGS = -g -O2 -fPIC/g' > temp
mv temp  ./src/packages/ampl/Makefile 

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

#Weird behaviour of Dakota: libamplsolver.a and amplsolver.a are not the same thing!
cd install/lib
mv libamplsolver.a libamplsolver.a.bak
ln -s ../../src/packages/ampl/amplsolver.a ./libamplsolver.a
