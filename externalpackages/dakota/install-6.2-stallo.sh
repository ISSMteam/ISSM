#!/bin/bash
set -eu


## Constants
#
VER="6.2"

#Some cleanup
rm -rf Dakota
rm -rf src 
rm -rf build 
rm -rf install 
mkdir src build install 

#Download from ISSM server
#$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/dakota-6.2-public.src.tar.gz' 'dakota-6.2-public-src.tar.gz'

#Untar 
tar -zxvf dakota-6.2-public-src.tar.gz

#Move Dakota to src directory
mv dakota-6.2.0.src/* src
rm -rf dakota-6.2.0.src

#Set up Dakota cmake variables and config
#export PATH="/usr/bin":$PATH
export DAK_SRC=$ISSM_DIR/externalpackages/dakota/src
export DAK_BUILD=$ISSM_DIR/externalpackages/dakota/build
export MPIHOME=/global/hds/software/cpu/eb3/impi/5.0.3.048-iccifort-2015.3.187-GNU-4.9.3-2.25/

cp ${DAK_SRC}/cmake/BuildDakotaTemplate.cmake $DAK_SRC/cmake/BuildDakotaCustom.cmake

# Patch source
patch ${DAK_SRC}/cmake/BuildDakotaCustom.cmake configs/${VER}/BuildDakotaCustom.cmake.stallo.patch
patch ${DAK_SRC}/cmake/DakotaDev.cmake configs/${VER}/DakotaDev.cmake.patch
patch ${DAK_SRC}/CMakeLists.txt configs/${VER}/CMakeLists.txt.stallo.patch
patch ${DAK_SRC}/src/dakota_data_io.hpp configs/${VER}/src/dakota_data_io.hpp.patch
patch ${DAK_SRC}/src/NonDSampling.cpp configs/${VER}/NonDSampling.cpp.patch
patch ${DAK_SRC}/src/NonDLocalReliability.cpp configs/${VER}/NonDLocalReliability.cpp.patch
patch ${DAK_SRC}/packages/pecos/src/pecos_global_defs.hpp configs/${VER}/pecos_global_defs.hpp.patch
patch ${DAK_SRC}/packages/surfpack/src/surfaces/nkm/NKM_KrigingModel.cpp configs/${VER}/NKM_KrigingModel.patch
patch ${DAK_SRC}/packages/DDACE/src/Analyzer/MainEffectsExcelOutput.cpp configs/${VER}/MainEffectsExcelOutput.patch
patch ${DAK_SRC}/src/DakotaInterface.cpp configs/${VER}/DakotaInterface.patch

#Configure dakota
cd $DAK_BUILD

cmake -D CMAKE_C_COMPILER=/global/hds/software/cpu/eb3/impi/5.0.3.048-iccifort-2015.3.187-GNU-4.9.3-2.25/bin64/mpicc \
	   -D CMAKE_CXX_COMPILER=/global/hds/software/cpu/eb3/impi/5.0.3.048-iccifort-2015.3.187-GNU-4.9.3-2.25/bin64/mpicxx \
	   -D CMAKE_Fortran_COMPILER=gfortran \
		-DHAVE_ACRO=off \
		-DHAVE_JEGA=off \
		-C $DAK_SRC/cmake/BuildDakotaCustom.cmake \
		-C $DAK_SRC/cmake/DakotaDev.cmake \
		$DAK_SRC
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
