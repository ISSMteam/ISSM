# Modifies path-related environment variables based on which external packages
# have been installed.
#
# ISSM_DIR and ISSM_ARCH should have been defined already in your shell
# settings file (i.e. .bashrc, .cshrc).
#
# TODO:
# - Condition all path modifications on existence of external package 'install'
#	directory
#

if [[ -n "$ZSH_VERSION" ]]; then
	# Silence `zsh: no matches found: <file>`
	setopt +o nomatch 1> /dev/null 2>&1
fi

## Functions
#
c_include_path_append(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${C_INCLUDE_PATH}" ]; then
			export C_INCLUDE_PATH="${1}"
		elif [[ ":${C_INCLUDE_PATH}:" != *":${1}:"* ]]; then
			export C_INCLUDE_PATH="${C_INCLUDE_PATH}:${1}"
		fi
	fi
} #}}}
c_include_path_prepend(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${C_INCLUDE_PATH}" ]; then
			export C_INCLUDE_PATH="${1}"
		elif [[ ":${C_INCLUDE_PATH}:" != *":${1}:"* ]]; then
			export C_INCLUDE_PATH="${1}:${C_INCLUDE_PATH}"
		fi
	fi
} #}}}

cpath_append(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${CPATH}" ]; then
			export CPATH="${1}"
		elif [[ ":${CPATH}:" != *":${1}:"* ]]; then
			export CPATH="${CPATH}:${1}"
		fi
	fi
} #}}}
cpath_prepend(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${CPATH}" ]; then
			export CPATH="${1}"
		elif [[ ":${CPATH}:" != *":${1}:"* ]]; then
			export CPATH="${1}:${CPATH}"
		fi
	fi
} #}}}

cplus_include_path_append(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${CPLUS_INCLUDE_PATH}" ]; then
			export CPLUS_INCLUDE_PATH="${1}"
		elif [[ ":${CPLUS_INCLUDE_PATH}:" != *":${1}:"* ]]; then
			export CPLUS_INCLUDE_PATH="${CPLUS_INCLUDE_PATH}:${1}"
		fi
	fi
} #}}}
cplus_include_path_prepend(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${CPLUS_INCLUDE_PATH}" ]; then
			export CPLUS_INCLUDE_PATH="${1}"
		elif [[ ":${CPLUS_INCLUDE_PATH}:" != *":${1}:"* ]]; then
			export CPLUS_INCLUDE_PATH="${1}:${CPLUS_INCLUDE_PATH}"
		fi
	fi
} #}}}

dyld_library_path_append(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${DYLD_LIBRARY_PATH}" ]; then
			export DYLD_LIBRARY_PATH="${1}"
		elif [[ ":${DYLD_LIBRARY_PATH}:" != *":${1}:"* ]]; then
			export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${1}"
		fi
		if [ -z "${LD_RUN_PATH}" ]; then
			export LD_RUN_PATH=$1
		elif [[ ":${LD_RUN_PATH}:" != *":${1}:"* ]]; then
			export LD_RUN_PATH="${LD_RUN_PATH}:${1}"
		fi
	fi
} #}}}
dyld_library_path_prepend(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${DYLD_LIBRARY_PATH}" ]; then
			export DYLD_LIBRARY_PATH="${1}"
		elif [[ ":${DYLD_LIBRARY_PATH}:" != *":${1}:"* ]]; then
			export DYLD_LIBRARY_PATH="${1}:${DYLD_LIBRARY_PATH}"
		fi
		if [ -z "${LD_RUN_PATH}" ]; then
			export LD_RUN_PATH="${1}"
		elif [[ ":${LD_RUN_PATH}:" != *":${1}:"* ]]; then
			export LD_RUN_PATH="${1}:${LD_RUN_PATH}"
		fi
	fi
} #}}}

ld_library_path_append(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${LD_LIBRARY_PATH}" ]; then
			export LD_LIBRARY_PATH="${1}"
		elif [[ ":${LD_LIBRARY_PATH}:" != *":${1}:"* ]]; then
			export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${1}"
		fi
		if [ -z "${LD_RUN_PATH}" ]; then
			export LD_RUN_PATH="${1}"
		elif [[ ":${LD_RUN_PATH}:" != *":$1:"* ]]; then
			export LD_RUN_PATH="${LD_RUN_PATH}:${1}"
		fi
	fi
} #}}}
ld_library_path_prepend(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${LD_LIBRARY_PATH}" ]; then
			export LD_LIBRARY_PATH="${1}"
		elif [[ ":${LD_LIBRARY_PATH}:" != *":${1}:"* ]]; then
			export LD_LIBRARY_PATH="${1}:${LD_LIBRARY_PATH}"
		fi
		if [ -z "${LD_RUN_PATH}" ]; then
			export LD_RUN_PATH="${1}"
		elif [[ ":${LD_RUN_PATH}:" != *":${1}:"* ]]; then
			export LD_RUN_PATH="${1}:${LD_RUN_PATH}"
		fi
	fi
} #}}}

library_path_append(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${LIBRARY_PATH}" ]; then
			export LIBRARY_PATH="${1}"
		elif [[ ":${LIBRARY_PATH}:" != *":$1:"* ]]; then
			export LIBRARY_PATH="${LIBRARY_PATH}:${1}"
		fi
	fi
} #}}}
library_path_prepend(){ #{{{
	if [ -d "${1}" ]; then
		if [ -z "${LIBRARY_PATH}" ]; then
			export LIBRARY_PATH="${1}"
		elif [[ ":${LIBRARY_PATH}:" != *":$1:"* ]]; then
			export LIBRARY_PATH="${1}:${LIBRARY_PATH}"
		fi
	fi
} #}}}

path_append(){ #{{{
	if [ -d "${1}" ] && [[ ":${PATH}:" != *":${1}:"* ]]; then
		PATH_IN="${1}"
		if [[ "${ISSM_ARCH}" == "cygwin-intel" ]]; then
			PATH_IN=`cygpath -u "${1}"`
		fi
		export PATH="${PATH}:${PATH_IN}"
	fi
} #}}}
path_prepend(){ #{{{
	if [ -d "${1}" ] && [[ ":${PATH}:" != *":${1}:"* ]]; then
		PATH_IN="${1}"
		if [[ "${ISSM_ARCH}" == "cygwin-intel" ]]; then
			PATH_IN=`cygpath -u "${1}"`
		fi
		export PATH="${PATH_IN}:${PATH}"
	fi
} #}}}

path_append "${ISSM_DIR}/aux-config"
path_append "${ISSM_DIR}/scripts"

# Default path to external packages. Redefine this constant if they are 
# installed to a different directory. Alternatively, export it on the command 
# line or in a profile initialization file (that is why we check here if it is 
# set already).
#
if [ -z "${ISSM_EXT_DIR+x}" ]; then
	export ISSM_EXT_DIR="${ISSM_DIR}/externalpackages"
fi

#######################
# OS-specific handling
#######################

OS_NAME=$(uname -s)

## macOS
#
if [[ ${OS_NAME} == "Darwin" ]]; then
	BUILD_TOOL_PATH=$(xcode-select -p)

	if [[ ${BUILD_TOOL_PATH} == "/Applications/Xcode.app/Contents/Developer" ]]; then
		BUILD_TOOL_VER=$(xcodebuild -version | /usr/bin/grep "Xcode" | sed -e 's/Xcode //' | cut -d. -f1)
	elif [[ ${BUILD_TOOL_PATH} == "/Library/Developer/CommandLineTools" ]]; then
		BUILD_TOOL_VER=$(pkgutil --pkg-info=com.apple.pkg.CLTools_Executables | /usr/bin/grep "version" | sed -e 's/version: //' | cut -d. -f1)
	else
		echo "Error: On macOS, either Xcode or the Command Line Tools must be installed!"
		exit 1
	fi

	if [[ ${BUILD_TOOL_VER} -ge 15 ]]; then
		 # Add flag for classic linker only if it has not already been added
		if [[ -z $(echo $LDFLAGS | /usr/bin/grep "\-Wl,\-ld_classic") ]]; then
			export LDFLAGS="${LDFLAGS} -Wl,-ld_classic"
		fi
	else
		export LDFLAGS="${LDFLAGS}" # At least set LDFLAGS to null string if it is not already set (used in installation of ISSM and some external packages). Can remove this when issues with new Xcode linker have settled.
	fi
fi

## Windows
#
MINGW=0
if [[ ${OS_NAME} == MINGW* ]]; then
	MINGW=1
	MSMPI_ROOT="${ISSM_EXT_DIR}/msmpi/install"
	if [ -d "${MSMPI_ROOT}" ]; then
		export MSMPI_ROOT # Used in installation of ParMETIS, ScaLAPACK
		cpath_prepend "${MSMPI_ROOT}/include"
		library_path_prepend "${MSMPI_ROOT}/lib"
	fi

	MPIEXEC_DIR=$(cygpath -u $(cygpath -ms "/c/Program Files/Microsoft MPI/Bin"))
	if [ -d "${MPIEXEC_DIR}" ]; then
		export MPIEXEC_DIR
		path_append "${MPIEXEC_DIR}"
	fi

	path_prepend "${ISSM_DIR}/bin" # Allows dynamic loader to find DLLs
fi

#########################
# Command-line utilities
#########################
SSH_ROOT="${ISSM_EXT_DIR}/ssh"
if [ -d "${SSH_ROOT}" ]; then
	path_append "${SSH_ROOT}"
fi

SVN_ROOT="${ISSM_EXT_DIR}/svn/install"
if [ -d "${SVN_ROOT}" ]; then
	path_prepend "${SVN_ROOT}/bin"
	ld_library_path_append "${SVN_ROOT}/lib"
fi

GIT_ROOT="${ISSM_EXT_DIR}/git/install"
if [ -d "${GIT_ROOT}" ]; then
	path_prepend "${GIT_ROOT}/bin"
fi

#############################
# Compilers / runtime / SDKs
#############################
export ANDROID_DIR="${ISSM_EXT_DIR}/android"

export ANDROID_NDK_DIR="$ANDROID_DIR/android-ndk/install"
path_append "$ANDROID_NDK_DIR/arm-linux-android-install/bin"

export ANDROID_SDK_ROOT="$ANDROID_DIR/android-sdk/install"
path_append "$ANDROID_SDK_ROOT/"

EMSCRIPTEN_ROOT="${ISSM_EXT_DIR}/emscripten/install"
if [ -d ${EMSCRIPTEN_ROOT} ]; then
	export EMSCRIPTEN_ROOT # Used in JavaScript build in installation of GSL, Triangle
fi

JVM_ROOT="/usr/local/gcc/4.3.2/lib64/gcj-4.3.2-9/"
ld_library_path_append "${JVM_ROOT}"

#############################
# IDEs
#############################
ECLIPSE_ROOT="${ISSM_EXT_DIR}/eclipse/install"
path_append "${ECLIPSE_ROOT}"

#############################
# Build systems
#############################
AUTOTOOLS_ROOT="${ISSM_EXT_DIR}/autotools/install"
path_prepend "${AUTOTOOLS_ROOT}/bin"

CMAKE_ROOT="${ISSM_EXT_DIR}/cmake/install"
path_prepend "${CMAKE_ROOT}/bin"

GMAKE_ROOT="${ISSM_EXT_DIR}/gmake/install"
path_prepend "${GMAKE_ROOT}/bin"

#############################
# Packagers
#############################
PACKAGEMAKER_ROOT="${ISSM_EXT_DIR}/packagemaker/install"
path_append "${PACKAGEMAKER_ROOT}"

#############################
# Libraries / binaries
#############################
MPI_ROOT_TEMP="${ISSM_EXT_DIR}/mpich/install"
if [ -d "${MPI_ROOT_TEMP}" ]; then
	export MPI_ROOT=${MPI_ROOT_TEMP}
	export MPI_DIR=${MPI_ROOT_TEMP}
	export MPI_HOME=${MPI_ROOT_TEMP} # Used in installation of Dakota
	export MPI_INC_DIR="${MPI_ROOT_TEMP}/include"
	path_prepend "${MPI_ROOT_TEMP}/bin"
	cpath_prepend "${MPI_ROOT_TEMP}/include"
	ld_library_path_append "${MPI_ROOT_TEMP}/lib"
fi

# NOTE: The following checks *must* come before PETSc as we prefer packages 
#		installed via PETSc
#
BLAS_ROOT="${ISSM_EXT_DIR}/blas/install"
if [ -d "${BLAS_ROOT}" ]; then
	export BLAS_ROOT # Used in installation of LAPACK, ScaLAPACK, PETSc
	library_path_append "${BLAS_ROOT}/lib"
	ld_library_path_append "${BLAS_ROOT}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${BLAS_ROOT}/lib" # Allows dynamic loader to find DLLs
	fi
fi

HDF5_ROOT="${ISSM_EXT_DIR}/hdf5/install"
if [ -d "${HDF5_ROOT}" ]; then
	export HDF5_ROOT # Used in installation of NetCDF, GDAL
	cpath_append "${HDF5_ROOT}/include"
	library_path_prepend "${HDF5_ROOT}/lib"
	dyld_library_path_prepend "${HDF5_ROOT}/lib"
	ld_library_path_prepend "${HDF5_ROOT}/lib"
fi

LAPACK_ROOT="${ISSM_EXT_DIR}/lapack/install"
if [ -d "${LAPACK_ROOT}" ]; then
	export LAPACK_ROOT # Used in installation of ScaLAPACK, MUMPS, PETSc
	library_path_append "${LAPACK_ROOT}/lib"
	ld_library_path_append "${LAPACK_ROOT}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${LAPACK_ROOT}/lib" # Allows dynamic loader to find DLLs
	fi

	if ls ${LAPACK_ROOT}/lib/libblas.* 1> /dev/null 2>&1; then
		export BLAS_ROOT="${LAPACK_ROOT}"
	fi
fi

METIS_ROOT="${ISSM_EXT_DIR}/metis/install"
if [ -d "${METIS_ROOT}" ]; then
	export METIS_ROOT # Used in installation of ParMETIS, Gmsh, PETSc
	library_path_prepend "${METIS_ROOT}/lib"
	ld_library_path_prepend "${METIS_ROOT}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${METIS_ROOT}/lib" # Allows dynamic loader to find DLLs
	fi
fi

MUMPS_ROOT="${ISSM_EXT_DIR}/mumps/install"
if [ -d "${MUMPS_ROOT}" ]; then
	export MUMPS_ROOT # Used in installation of PETSc
	library_path_append "${MUMPS_ROOT}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${MUMPS_ROOT}/lib" # Allows dynamic loader to find DLLs
	fi
fi

PARMETIS_ROOT="${ISSM_EXT_DIR}/parmetis/install"
if [ -d "${PARMETIS_ROOT}" ]; then
	export PARMETIS_ROOT # Used in installation of MUMPS, PETSc
	library_path_prepend "${PARMETIS_ROOT}/lib"
	ld_library_path_prepend "${PARMETIS_ROOT}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${PARMETIS_ROOT}/lib" # Allows dynamic loader to find DLLs
	fi
fi

QD_ROOT="${ISSM_EXT_DIR}/qd/install"
if [ -d "${QD_ROOT}" ]; then
	export QD_ROOT # Used in installation of MPLAPACK
	library_path_prepend "${QD_ROOT}/lib"
	ld_library_path_prepend "${QD_ROOT}/lib"
fi

SCALAPACK_ROOT="${ISSM_EXT_DIR}/scalapack/install"
if [ -d "${SCALAPACK_ROOT}" ]; then
	export SCALAPACK_ROOT # Used in installation of MUMPS, PETSc
	library_path_append "${SCALAPACK_ROOT}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${SCALAPACK_ROOT}/lib" # Allows dynamic loader to find DLLs
	fi
fi

ZLIB_ROOT="${ISSM_EXT_DIR}/zlib/install"
if [ -d "${ZLIB_ROOT}" ]; then
	export ZLIB_ROOT # Used in installation of NetCDF, GDAL, GMT
	ld_library_path_append "${ZLIB_ROOT}/lib"
fi

PETSC_ROOT="${ISSM_EXT_DIR}/petsc/install"
if [ -d "${PETSC_ROOT}" ]; then
	export PETSC_ROOT # Used in installation of Gmsh
	cpath_prepend "${PETSC_ROOT}/include"
	library_path_prepend "${PETSC_ROOT}/lib"
	ld_library_path_prepend "${PETSC_ROOT}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${PETSC_ROOT}/lib" # Allows dynamic loader to find DLLs
	fi

	# In case we have installed certain external packages via PETSc
	#

	# BLAS
	if ls ${PETSC_ROOT}/lib/libblas.* 1> /dev/null 2>&1 || ls ${PETSC_ROOT}/lib/libfblas.* 1> /dev/null 2>&1; then
		export BLAS_ROOT="${PETSC_ROOT}" # Used in installation of Dakota, GMT
	fi

	# HDF5
	if ls ${PETSC_ROOT}/lib/libhdf5.* 1> /dev/null 2>&1; then
		export HDF5_ROOT="${PETSC_ROOT}" # Used in installation of NetCDF, GDAL
		cpath_append "${PETSC_ROOT}/include"
		library_path_append "${PETSC_ROOT}/lib"
		dyld_library_path_append "${PETSC_ROOT}/lib"
		ld_library_path_append "${PETSC_ROOT}/lib"
	fi

	# LAPACK
	if ls ${PETSC_ROOT}/lib/liblapack.* 1> /dev/null 2>&1 || ls ${PETSC_ROOT}/lib/libflapack.* 1> /dev/null 2>&1; then
		export LAPACK_ROOT="${PETSC_ROOT}" # Used in installation of Dakota, GMT
	fi

	# METIS
	if ls ${PETSC_ROOT}/lib/libmetis.* 1> /dev/null 2>&1; then
		export METIS_ROOT="${PETSC_ROOT}" # Used in installation of Gmsh
	fi

	# MPICH
	if [ -f "${PETSC_ROOT}/bin/mpiexec" ]; then
		export MPI_ROOT=${PETSC_ROOT}
		export MPI_DIR=${MPI_ROOT}
		export MPI_HOME=${MPI_ROOT} # Used in installation of Dakota
		export MPI_INC_DIR="${MPI_ROOT}/include"
		path_prepend "${MPI_ROOT}/bin"
		cpath_prepend "${MPI_ROOT}/include"
		ld_library_path_append "${MPI_ROOT}/lib"
	fi

	# ZLIB
	if ls ${PETSC_ROOT}/lib/libz.* 1> /dev/null 2>&1; then
		export ZLIB_ROOT="${PETSC_ROOT}" # Used in installation of NetCDF, GDAL
		ld_library_path_append "${PETSC_ROOT}/lib"
	fi
fi

MPLAPACK_ROOT="${ISSM_EXT_DIR}/mplapack/install"
if [ -d "${MPLAPACK_ROOT}" ]; then
	cplus_include_path_prepend "${MPLAPACK_ROOT}/include"
	cplus_include_path_prepend "${MPLAPACK_ROOT}/include/mplapack"
	library_path_prepend "${MPLAPACK_ROOT}/lib"
	ld_library_path_prepend "${MPLAPACK_ROOT}/lib"
fi

SCOTCH_ROOT="${ISSM_EXT_DIR}/scotch/install"
ld_library_path_append "${SCOTCH_ROOT}/lib"

SLEPC_ROOT="${ISSM_EXT_DIR}/slepc/install"
ld_library_path_append "${SLEPC_ROOT}/lib"

TAO_ROOT="${ISSM_EXT_DIR}/tao/install"
ld_library_path_append "${TAO_ROOT}/lib"

BOOST_ROOT="${ISSM_EXT_DIR}/boost/install"
if [ -d "${BOOST_ROOT}" ]; then
	export BOOST_ROOT # Used in installation of Dakota
	export BOOST_DIR=${BOOST_ROOT}
	export BOOSTROOT=${BOOST_ROOT}
	path_append "${BOOST_ROOT}/bin"
	library_path_prepend "${BOOST_ROOT}/lib"
	ld_library_path_prepend "${BOOST_ROOT}/lib"
	dyld_library_path_prepend "${BOOST_ROOT}/lib"
fi

DAKOTA_ROOT="${ISSM_EXT_DIR}/dakota/install"
if [ -d "${DAKOTA_ROOT}" ]; then
	path_append "${DAKOTA_ROOT}/bin"
	ld_library_path_prepend "${DAKOTA_ROOT}/lib"
	dyld_library_path_prepend "${DAKOTA_ROOT}/lib"
fi

NCO_ROOT="${ISSM_EXT_DIR}/nco/install/bin"
path_prepend "${NCO_ROOT}/bin"

CPPCHECK_ROOT="${ISSM_EXT_DIR}/cppcheck/install"
path_append "${CPPCHECK_ROOT}/bin"

MERCURIAL_ROOT="${ISSM_EXT_DIR}/mercurial/install"
if [ -d "${MERCURIAL_ROOT}" ]; then
	export PYTHONPATH="${PYTHONPATH}:${MERCURIAL_ROOT}/mercurial/pure/"
	path_append "${MERCURIAL_ROOT}"
fi

GSL_ROOT="${ISSM_EXT_DIR}/gsl/install"
if [ -d "${GSL_ROOT}" ]; then
	export GSL_HOME="${GSL_ROOT}" # Used in installation of Dakota
	cpath_prepend "${GSL_ROOT}/include"
	ld_library_path_append "${GSL_ROOT}/lib"
fi

NETCDF_ROOT="${ISSM_EXT_DIR}/netcdf/install"
if [ -d "${NETCDF_ROOT}" ]; then
	export NETCDF_ROOT # Used in installation of GDAL, GMT
	path_prepend "${NETCDF_ROOT}/bin"
	cpath_prepend "${NETCDF_ROOT}/include"
	library_path_prepend "${NETCDF_ROOT}/lib"
	ld_library_path_prepend "${NETCDF_ROOT}/lib"
	dyld_library_path_prepend "${NETCDF_ROOT}/lib"
fi

NETCDF_CXX_ROOT="${ISSM_EXT_DIR}/netcdf-cxx/install"
if [ -d "${NETCDF_CXX_ROOT}" ]; then
	ld_library_path_append "${NETCDF_CXX_ROOT}/lib"
fi

NETCDF_PYTHON_ROOT="${ISSM_EXT_DIR}/netcdf-python/install"
if [ -d "${NETCDF_PYTHON_ROOT}" ]; then
	if [ -d "${NETCDF_PYTHON_ROOT}/lib/python2.7/site-packages" ]; then
		ld_library_path_append "${NETCDF_PYTHON_ROOT}/lib/python2.7/site-packages"
	fi
fi

CURL_ROOT="${ISSM_EXT_DIR}/curl/install"
if [ -d "${CURL_ROOT}" ]; then
	export CURL_ROOT # Used in installation of NetCDF, GDAL, GMT
	cpath_prepend "${CURL_ROOT}/include"
	ld_library_path_prepend "${CURL_ROOT}/lib"
	dyld_library_path_prepend "${CURL_ROOT}/lib"
	path_append "${CURL_ROOT}/bin"
fi

SQLITE_ROOT="${ISSM_EXT_DIR}/sqlite/install"
if [ -d "${SQLITE_ROOT}" ]; then
	export SQLITE_ROOT # Used in installation of GDAL
	path_prepend "${SQLITE_ROOT}/bin"
	cpath_prepend "${SQLITE_ROOT}/include"
	library_path_prepend "${SQLITE_ROOT}/lib"
	ld_library_path_prepend "${SQLITE_ROOT}/lib"
fi

LIBTIFF_ROOT="${ISSM_EXT_DIR}/libtiff/install"
if [ -d "${LIBTIFF_ROOT}" ]; then
	dyld_library_path_append "${LIBTIFF_ROOT}/install/libtiff"
	ld_library_path_append "${LIBTIFF_ROOT}/install/libtiff"
fi

PROJ_ROOT="${ISSM_EXT_DIR}/proj/install"
if [ -d "${PROJ_ROOT}" ]; then
	export PROJ_ROOT # Used in installation of GDAL
	path_append "${PROJ_ROOT}/bin"
	ld_library_path_prepend "${PROJ_ROOT}/lib"
	dyld_library_path_prepend "${PROJ_ROOT}/lib"
	ld_library_path_prepend "${PROJ_ROOT}/lib"
fi

GDAL_ROOT="${ISSM_EXT_DIR}/gdal/install"
if [ -d "${GDAL_ROOT}" ]; then
	export GDAL_ROOT # Used in installation of GMT
	path_prepend "${GDAL_ROOT}/bin"
	ld_library_path_append "${GDAL_ROOT}/lib"
	dyld_library_path_append "${GDAL_ROOT}/lib"
fi

GSHHG_ROOT="${ISSM_EXT_DIR}/gshhg/install"
if [ -d "${GSHHG_ROOT}" ]; then
	export GSHHG_ROOT # Used in installation of GMT
fi

GMT_ROOT="${ISSM_EXT_DIR}/gmt/install"
if [ -d "${GMT_ROOT}" ]; then
	path_prepend "${GMT_ROOT}/bin"
	ld_library_path_append "${GMT_ROOT}/lib"
	dyld_library_path_append "${GMT_ROOT}/lib"
fi

GMSH_ROOT="${ISSM_EXT_DIR}/gmsh/install"
if [ -d "${GMSH_ROOT}" ]; then
	path_prepend "${GMSH_ROOT}/bin"
	ld_library_path_append "${GMSH_ROOT}/lib"
	dyld_library_path_append "${GMSH_ROOT}/lib"
fi

TRIANGLE_ROOT="${ISSM_EXT_DIR}/triangle/install"
if [ -d "${TRIANGLE_ROOT}" ]; then
	ld_library_path_append "${TRIANGLE_ROOT}/lib"
	dyld_library_path_append "${TRIANGLE_ROOT}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${TRIANGLE_ROOT}/lib" # Allows dynamic loader to find DLLs
	fi
fi

ANGELROOT="${ISSM_EXT_DIR}/angel/angel"
if [ -d "${ANGELROOT}" ]; then
	export ANGELROOT
fi

OPENANALYSISROOT="${ISSM_EXT_DIR}/openanalysis/install"
if [ -d "${OPENANALYSISROOT}" ]; then
	export OPENANALYSISROOT
	ld_library_path_append "${OPENANALYSISROOT}/lib"
fi

BBFTP_ROOT="${ISSM_EXT_DIR}/bbftp/install"
path_append "${BBFTP_ROOT}/bin"

ADIC_ROOT="${ISSM_EXT_DIR}/adic/install"
path_append "${ADIC_ROOT}/bin"
ld_library_path_append "${ADIC_ROOT}/lib"

COLPACK_ROOT="${ISSM_EXT_DIR}/colpack/install"
ld_library_path_append "${COLPACK_ROOT}/lib"

APPSCAN_ROOT="${ISSM_EXT_DIR}/appscan/install"
path_append "${APPSCAN_ROOT}/bin"

RATS_ROOT="${ISSM_EXT_DIR}/rats/install"
path_append "${RATS_ROOT}/bin"

DYSON_ROOT="${ISSM_EXT_DIR}/dyson/"
path_append "${DYSON_ROOT}"

SHAPELIB_ROOT="${ISSM_EXT_DIR}/shapelib/install"
path_append "${SHAPELIB_ROOT}/exec"

CCCL_ROOT="${ISSM_EXT_DIR}/cccl/install"
path_append "${CCCL_ROOT}/bin"

MODELE_ROOT="${ISSM_EXT_DIR}/modelE/install"
path_append "${MODELE_ROOT}/src/exec"

NCVIEW_ROOT="${ISSM_EXT_DIR}/ncview/install"
path_append "${NCVIEW_ROOT}"

TCLX_ROOT="${ISSM_EXT_DIR}/tclx/install/lib/tclx8.4"
ld_library_path_append "${TCLX_ROOT}"

ASPELL_ROOT="${ISSM_EXT_DIR}/aspell/install"
path_append "${ASPELL_ROOT}/bin"

ESMF_ROOT="${ISSM_EXT_DIR}/esmf/install"
if [ -d "${ESMF_ROOT}" ]; then
	path_prepend "${ESMF_ROOT}/bin"
	ld_library_path_append "${ESMF_ROOT}/lib/libO/Linux.gfortran.64.mpich.default/"
fi

CVS_ROOT="${ISSM_EXT_DIR}/cvs/install"
path_prepend "${CVS_ROOT}/bin"

APR_ROOT="${ISSM_EXT_DIR}/apr/install"
path_append "${APR_ROOT}/bin"
ld_library_path_append "${APR_ROOT}/lib"

APR_UTIL_ROOT="${ISSM_EXT_DIR}/apr-util/install"
path_prepend "${APR_UTIL_ROOT}/bin"
ld_library_path_append "${APR_UTIL_ROOT}/lib"

YAMS_ROOT="${ISSM_EXT_DIR}/yams/install"
path_append "${YAMS_ROOT}"

SWIG_ROOT="${ISSM_EXT_DIR}/swig/install"
path_append "${SWIG_ROOT}"

INISHELL_ROOT="${ISSM_EXT_DIR}/inishell/install"
path_append "${INISHELL_ROOT}"

EXPAT_ROOT="${ISSM_EXT_DIR}/expat/install"
ld_library_path_prepend "${EXPAT_ROOT}"
dyld_library_path_prepend "${EXPAT_ROOT}"

NEOPZ_ROOT="${ISSM_EXT_DIR}/neopz/install"
if [ -d "${NEOPZ_ROOT}" ]; then
	export REFPATTERNDIR="${NEOPZ_ROOT}/include/refpatterns"
fi

XERCESROOT="${ISSM_EXT_DIR}/xerces/install"
if [ -d "${XERCESROOT}" ]; then
	export XERCESROOT
	export XERCESCROOT="${ISSM_EXT_DIR}/xerces/src"
fi

XAIFBOOSTERROOT="${ISSM_EXT_DIR}/xaifbooster"
XAIF_ROOT="${XAIFBOOSTERROOT}/xaifBooster"
if [ -d "${XAIF_ROOT}" ]; then
	export XAIFBOOSTERROOT
	export XAIF_DIR="${XAIF_ROOT}"
	export XAIFBOOSTER_HOME="${XAIF_ROOT}"
	export PLATFORM="x86-Linux"
fi

VALGRIND_ROOT="${ISSM_DIR}/externalpackages/valgrind/install"
if [ -d "${VALGRIND_ROOT="${ISSM_DIR}/valgrind/install"
}" ]; then
	path_prepend "${VALGRIND_ROOT}/bin"
fi

DOXYGEN_ROOT="${ISSM_EXT_DIR}/doxygen/install"
path_prepend "${DOXYGEN_ROOT}/bin"

SHELL2JUNIT_ROOT="${ISSM_EXT_DIR}/shell2junit/install"
path_append "${SHELL2JUNIT_ROOT}"
