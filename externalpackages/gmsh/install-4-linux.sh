#!/bin/bash
set -e


# TODO:
# - Add support for,
#	- MUMPS
#	- PETSc
# 
# See Also:
# - configs/4/static/CMakeLists.txt
# - http://gmsh.info/doc/texinfo/gmsh.html#Compiling-the-source-code
#
# NOTE: Currently downloading master development branch to overcome the issue 
# described here: https://gitlab.onelab.info/gmsh/gmsh/-/issues/3158. Revert 
# changes once version 4.13.2 has been released.
#

## Constants
#
VER="4.12.2"

PREFIX="${ISSM_DIR}/externalpackages/gmsh/install" # Set to location where external package should be installed

# Find libgfortran so that we do not have to hardcode it.
#
# TODO:
# - Move this to etc/environment.sh
#
echo "Finding libgfortran..."
LIBGFORTRAN=$(find /usr -name libgfortran.* 2>/dev/null | egrep -n libgfortran.a | sed "s/[0-9]*://g" | head -1)
LIBGFORTRAN_ROOT=${LIBGFORTRAN%/*}

# Environment
#
export CXXFLAGS="${CXXFLAGS} -w"

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} src

# Download source
#$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://gmsh.info/src/gmsh-${VER}-source.tgz" "gmsh-${VER}-source.tgz"
# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://gitlab.onelab.info/gmsh/gmsh/-/archive/master/gmsh-master.tar.gz" "gmsh-master.tar.gz"

# Unpack source
#tar -xvzf gmsh-${VER}-source.tgz
tar -xvzf gmsh-master.tar.gz

# Move source to 'src' directory
#mv gmsh-${VER}-source/* src
#rm -rf gmsh-${VER}-source
mv gmsh-master/* src
rm -rf gmsh-master

# Configure
#
# NOTE:
# - Option -DENABLE_FLTK=0 is used because we do not need GUI.
# - Option -DENABLE_MPEG_ENCODE=0 is used because we do not need to record MPEG 
#	movies.
# - Option -DENABLE_OCC=0 is used because we do not need CAD kernel and are not 
#	importing STEP/IGES files.
# - Option -DENABLE_TOUCHBAR=0 is used because we do not have GUI, therefore we 
#	do not need to support Apple Touch bar (does not affect compilation on 
#	Linux).
#
cd src
cmake \
	-DCMAKE_INSTALL_PREFIX="${PREFIX}" \
	-DCMAKE_INSTALL_RPATH="${PREFIX}/lib" \
	-DCMAKE_BUILD_TYPE=Release \
	-DENABLE_BUILD_DYNAMIC=1 \
	-DENABLE_BUILD_SHARED=1 \
	-DBLAS_LAPACK_LIBRARIES="-L${LAPACK_ROOT}/lib -lflapack -L${BLAS_ROOT}/lib -lfblas -L${LIBGFORTRAN_ROOT} -lgfortran" \
	-DENABLE_BLAS_LAPACK=1 \
	-DENABLE_EIGEN=0 \
	-DENABLE_FLTK=0 \
	-DENABLE_MPEG_ENCODE=0 \
	-DENABLE_MPI=1 \
	-DENABLE_OCC=0 \
	-DENABLE_TOUCHBAR=0 \
	-DMETIS_ROOT="${METIS_ROOT}"

# Compile and install
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi

# Make necessary link on RHEL
if [[ -d ${PREFIX}/lib64 && ! -d ${PREFIX}/lib ]]; then
	cd ${PREFIX}
	ln -s ./lib64 ./lib
fi
