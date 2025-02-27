#!/bin/bash
set -eu


## Constants
#
VER="6.2"

PREFIX="${ISSM_DIR}/externalpackages/dakota/install" # Set to location where external package should be installed

# Find libgfortran so that we do not have to hardcode it.
#
# TODO:
# - Move this to etc/environment.sh
#
echo "Finding libgfortran..."
LIBGFORTRAN=$(find /usr /opt -name libgfortran* 2>/dev/null | egrep -n libgfortran.a | egrep -v i386 | sed "s/[0-9]*://g" | head -1)
LIBGFORTRAN_ROOT=${LIBGFORTRAN%/*}

## Environment
#
export BLAS_LIBS="-L${BLAS_ROOT}/lib -lfblas" # Need to export BLAS_LIBS *and* pass it as an option to CMake to ensure that external packages also find it
export DAK_BUILD=${ISSM_DIR}/externalpackages/dakota/build # DO NOT CHANGE THIS
export DAK_INSTALL=${PREFIX} # DO NOT CHANGE THIS
export DAK_SRC=${ISSM_DIR}/externalpackages/dakota/src # DO NOT CHANGE THIS
export FLIBS="-L${LIBGFORTRAN_ROOT} -lgfortran"
export LAPACK_LIBS="-L${LAPACK_ROOT}/lib -lflapack" # Need to export LAPACK_LIBS *and* pass it as an option to CMake to ensure that external packages also find it
export LDFLAGS="-framework CoreFoundation"

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
cp configs/${VER}/cmake/BuildDakotaCustom.cmake ${DAK_SRC}/cmake
cp configs/${VER}/cmake/DakotaDev.cmake ${DAK_SRC}/cmake

# Copy customized source and configuration files specific to Mac to 'src' directory
cp configs/${VER}/mac/cmake/InstallDarwinDylibs.cmake ${DAK_SRC}/cmake

# Patch source
patch ${DAK_SRC}/src/dakota_data_io.hpp configs/${VER}/src/dakota_data_io.hpp.patch

# Disable requirement of Python 2 for TriBITS
sed -i'' -e 's|SET(PythonInterp_FIND_VERSION|#SET(PythonInterp_FIND_VERSION|' ${DAK_SRC}/packages/teuchos/cmake/tribits/package_arch/TribitsFindPythonInterp.cmake

# Configure
#
# NOTE:
# - The -w option has been added to CMAKE_C_FLAGS, CMAKE_CXX_FLAGS, and 
#	CMAKE_Fortran_FLAGS. This should be removed for more recent versions of 
#	Dakota.
#
cd ${DAK_BUILD}
cmake \
	-DBUILD_SHARED_LIBS=ON \
	-DBUILD_STATIC_LIBS=OFF \
	-DCMAKE_C_COMPILER=${MPI_HOME}/bin/mpicc \
	-DCMAKE_C_FLAGS="-w -Wno-error=implicit-int -Wno-implicit-function-declaration" \
	-DCMAKE_CXX_COMPILER=${MPI_HOME}/bin/mpicxx \
	-DCMAKE_CXX_FLAGS="-fdelayed-template-parsing -w" \
	-DCMAKE_CXX_STANDARD="11" \
	-DCMAKE_Fortran_COMPILER=${MPI_HOME}/bin/mpif77 \
	-DCMAKE_Fortran_FLAGS="-fallow-argument-mismatch -w" \
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

# Set install_name for all shared libraries
cd ${DAK_INSTALL}/lib
for name in *.dylib; do
	install_name_tool -id ${DAK_INSTALL}/lib/${name} ${name}
done

## Patch install names for certain libraries
#
# TODO: Figure out how to reconfigure source to apply these install names at 
# 		compile time
#
install_name_tool -change libdakota_src_fortran.dylib ${DAK_INSTALL}/lib/libdakota_src_fortran.dylib libdakota_src.dylib
install_name_tool -change liblhs_mod.dylib ${DAK_INSTALL}/lib/liblhs_mod.dylib liblhs.dylib
install_name_tool -change liblhs_mods.dylib ${DAK_INSTALL}/lib/liblhs_mods.dylib liblhs.dylib
install_name_tool -change liblhs_mod.dylib ${DAK_INSTALL}/lib/liblhs_mod.dylib liblhs_mods.dylib
install_name_tool -change libteuchos.dylib ${DAK_INSTALL}/lib/libteuchos.dylib liboptpp.dylib
install_name_tool -change libdfftpack.dylib ${DAK_INSTALL}/lib/libdfftpack.dylib libpecos.dylib
install_name_tool -change liblhs.dylib ${DAK_INSTALL}/lib/liblhs.dylib libpecos.dylib
install_name_tool -change liblhs_mod.dylib ${DAK_INSTALL}/lib/liblhs_mod.dylib libpecos.dylib
install_name_tool -change liblhs_mods.dylib ${DAK_INSTALL}/lib/liblhs_mods.dylib libpecos.dylib
install_name_tool -change libpecos_src.dylib ${DAK_INSTALL}/lib/libpecos_src.dylib libpecos.dylib
install_name_tool -change libteuchos.dylib ${DAK_INSTALL}/lib/libteuchos.dylib libpecos.dylib
install_name_tool -change libdfftpack.dylib ${DAK_INSTALL}/lib/libdfftpack.dylib libpecos_src.dylib
install_name_tool -change liblhs.dylib ${DAK_INSTALL}/lib/liblhs.dylib libpecos_src.dylib
install_name_tool -change liblhs_mod.dylib ${DAK_INSTALL}/lib/liblhs_mod.dylib libpecos_src.dylib
install_name_tool -change liblhs_mods.dylib ${DAK_INSTALL}/lib/liblhs_mods.dylib libpecos_src.dylib
install_name_tool -change libteuchos.dylib ${DAK_INSTALL}/lib/libteuchos.dylib libpecos_src.dylib
install_name_tool -change libsurfpack_fortran.dylib ${DAK_INSTALL}/lib/libsurfpack_fortran.dylib libsurfpack.dylib

## Add LIBGFORTRAN_ROOT to rpath for libraries that need it
#
# TODO: Figure out how to reconfigure source to add to rpath at compile time
#
install_name_tool -add_rpath ${LIBGFORTRAN_ROOT} libpecos.dylib
install_name_tool -add_rpath ${LIBGFORTRAN_ROOT} libteuchos.dylib
install_name_tool -add_rpath ${LIBGFORTRAN_ROOT} liboptpp.dylib

