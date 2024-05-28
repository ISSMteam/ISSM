#!/bin/bash
set -eu

#Some cleanup
rm -rf Dakota
rm -rf src 
rm -rf build 
rm -rf install 
mkdir src build install 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/dakota-5.3.1-public-src.tar.gz' 'dakota-5.3.1-public-src.tar.gz'

#Untar 
tar -zxvf dakota-5.3.1-public-src.tar.gz

#Move Dakota to src directory
mv dakota-5.3.1.src/* src
rm -rf dakota-5.3.1.src

#Set up Dakota cmake variables and config
export DAK_SRC=$ISSM_DIR/externalpackages/dakota/src
export DAK_BUILD=$ISSM_DIR/externalpackages/dakota/build
cp $DAK_SRC/cmake/BuildDakotaTemplate.cmake $DAK_SRC/cmake/BuildDakotaCustom.cmake
patch $DAK_SRC/cmake/BuildDakotaCustom.cmake configs/5.3.1/BuildDakotaCustom.cmake.mac.patch
patch $DAK_SRC/cmake/DakotaDev.cmake configs/5.3.1/DakotaDev.cmake.patch
patch $DAK_SRC/CMakeLists.txt configs/5.3.1/CMakeLists.txt.patch

#Apply patches
patch src/src/ParallelLibrary.cpp configs/5.3.1/ParallelLibrary.cpp.patch
patch src/src/ParallelLibrary.hpp configs/5.3.1/ParallelLibrary.hpp.patch
patch src/src/NonDSampling.cpp configs/5.3.1/NonDSampling.cpp.patch
patch src/src/NonDLocalReliability.cpp configs/5.3.1/NonDLocalReliability.cpp.patch
patch src/packages/pecos/src/pecos_global_defs.hpp configs/5.3.1/pecos_global_defs.hpp.patch

export BOOST_ROOT=$ISSM_DIR/externalpackages/boost/install

#Configure dakota
# Set your local gcc compiler here
cd $DAK_BUILD
cmake -DBoost_NO_BOOST_CMAKE=TRUE \
    -DBoost_NO_SYSTEM_PATHS=TRUE \
    -DBOOST_ROOT:PATHNAME=$BOOST_ROOT \
    -DBoost_LIBRARY_DIRS:FILEPATH=${BOOST_ROOT}/lib \
	 -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CC_COMPILER=gcc \
    -C $DAK_SRC/cmake/BuildDakotaCustom.cmake -C $DAK_SRC/cmake/DakotaDev.cmake $DAK_SRC

#-DCMAKE_CXX_COMPILER=/usr/local/gfortran/bin/x86_64-apple-darwin10-g++ -DCMAKE_Fortran_COMPILER=/usr/local/gfortran/bin/x86_64-apple-darwin10-gfortran

cd ..

#Compile and install dakota
cd $DAK_BUILD
if [ $# -eq 0 ];
then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
cd ..
