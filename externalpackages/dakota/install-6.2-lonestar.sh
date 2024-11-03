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
export DAK_SRC=$ISSM_DIR/externalpackages/dakota/src
export DAK_BUILD=$ISSM_DIR/externalpackages/dakota/build
export MPIHOME=/opt/cray/pe/mpt/7.7.3/gni/mpich-intel/16.0/

# Copy customized source and configuration files to 'src' directory
cp configs/${VER}/packages/DDACE/src/Analyzer/MainEffectsExcelOutput.cpp ${DAK_SRC}/packages/DDACE/src/Analyzer
cp configs/${VER}/packages/queso/src/misc/src/1DQuadrature.C ${DAK_SRC}/packages/queso/src/misc/src
cp configs/${VER}/packages/surfpack/src/surfaces/nkm/NKM_KrigingModel.cpp ${DAK_SRC}/packages/surfpack/src/surfaces/nkm
cp configs/${VER}/packages/VPISparseGrid/src/sandia_rules.cpp ${DAK_SRC}/packages/VPISparseGrid/src
cp configs/${VER}/src/DakotaInterface.cpp ${DAK_SRC}/src
cp configs/${VER}/src/NonDLocalReliability.cpp ${DAK_SRC}/src
cp configs/${VER}/src/NonDSampling.cpp ${DAK_SRC}/src
cp configs/${VER}/cmake/BuildDakotaCustom.cmake ${DAK_SRC}/cmake
cp configs/${VER}/cmake/DakotaDev.cmake ${DAK_SRC}/cmake

# Patch source
patch ${DAK_SRC}/src/dakota_data_io.hpp configs/${VER}/src/dakota_data_io.hpp.patch
patch ${DAK_SRC}/CMakeLists.txt configs/${VER}/CMakeLists.txt.lonestar.patch
patch ${DAK_SRC}/src/DakotaInterface.cpp configs/${VER}/DakotaInterface.patch



#Configure dakota
cd $DAK_BUILD

cmake -D CMAKE_C_COMPILER=mpicc \
	   -D CMAKE_CXX_COMPILER=mpicxx \
	   -D CMAKE_Fortran_COMPILER=mpif77 \
		-D MPIEXEC_EXECUTABLE=/opt/apps/tacc/bin/ibrun \
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
