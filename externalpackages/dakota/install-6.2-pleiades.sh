#!/bin/bash
set -eu


## Constants
#
VER="6.2"

PREFIX="${ISSM_DIR}/externalpackages/dakota/install" # Set to location where external package should be installed

## Environment
#
export BLAS_LIBS="-L${COMP_INTEL_ROOT}/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" # Need to export BLAS_LIBS *and* pass it as an option to CMake to ensure that external packages also find it
export CXXFLAGS='-std=c++98'
export DAK_BUILD=${ISSM_DIR}/externalpackages/dakota/build # DO NOT CHANGE THIS
export DAK_INSTALL=${PREFIX} # DO NOT CHANGE THIS
export DAK_SRC=${ISSM_DIR}/externalpackages/dakota/src # DO NOT CHANGE THIS
export LAPACK_LIBS="-L${COMP_INTEL_ROOT}/mkl/lib/intel64/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" # Need to export LAPACK_LIBS *and* pass it as an option to CMake to ensure that external packages also find it

# Cleanup
rm -rf ${DAK_BUILD} ${DAK_INSTALL} ${DAK_SRC}
mkdir -p ${DAK_BUILD} ${DAK_INSTALL} ${DAK_SRC}

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/dakota-${VER}-public.src.tar.gz" "dakota-${VER}-public-src.tar.gz"

# Unpack source
tar -zxvf dakota-${VER}-public-src.tar.gz

# Move source to 'src' directory
mv dakota-${VER}.0.src/* ${DAK_SRC}
rm -rf dakota-${VER}.0.src

# Copy customized source and configuration files to 'src' directory
cp configs/${VER}/packages/DDACE/src/Analyzer/MainEffectsExcelOutput.cpp ${DAK_SRC}/packages/DDACE/src/Analyzer
cp configs/${VER}/packages/surfpack/src/surfaces/nkm/NKM_KrigingModel.cpp ${DAK_SRC}/packages/surfpack/src/surfaces/nkm
cp configs/${VER}/packages/VPISparseGrid/src/sandia_rules.cpp ${DAK_SRC}/packages/VPISparseGrid/src
cp configs/${VER}/src/DakotaInterface.cpp ${DAK_SRC}/src
cp configs/${VER}/src/NonDLocalReliability.cpp ${DAK_SRC}/src
cp configs/${VER}/src/NonDSampling.cpp ${DAK_SRC}/src

# Copy customized source and configuration files specific to Linux to 'src' directory
cp configs/${VER}/linux/cmake/BuildDakotaCustom.pleiades.cmake ${DAK_SRC}/cmake/BuildDakotaCustom.cmake
cp configs/${VER}/cmake/DakotaDev.cmake ${DAK_SRC}/cmake

# Patch source
patch ${DAK_SRC}/src/dakota_data_io.hpp configs/${VER}/src/dakota_data_io.hpp.patch

# Disable requirement of Python 2 for TriBITS
sed -i'' -e 's|SET(PythonInterp_FIND_VERSION|#SET(PythonInterp_FIND_VERSION|' ${DAK_SRC}/packages/teuchos/cmake/tribits/package_arch/TribitsFindPythonInterp.cmake

# Configure
cd ${DAK_BUILD}
cmake \
	-DBUILD_SHARED_LIBS=ON \
	-DBUILD_STATIC_LIBS=OFF \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_C_FLAGS="-Wno-error=implicit-function-declaration" \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_Fortran_COMPILER=gfortran \
	-DBoost_NO_BOOST_CMAKE=TRUE \
	-DHAVE_ACRO=OFF \
	-DHAVE_JEGA=OFF \
	-C${DAK_SRC}/cmake/BuildDakotaCustom.cmake \
	-C${DAK_SRC}/cmake/DakotaDev.cmake \
	${DAK_SRC}

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi

cd ${DAK_INSTALL}

# Comment out definition of HAVE_MPI in Teuchos config header file in order to
# avoid conflict with our definition
sed -i -e "s/#define HAVE_MPI/\/* #define HAVE_MPI *\//g" include/Teuchos_config.h
