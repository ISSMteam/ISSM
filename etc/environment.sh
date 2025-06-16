# Modifies path-related environment variables based on which external packages
# have been installed.
#
# ISSM_DIR and ISSM_ARCH should have been defined already in your shell
# settings file (i.e., .bashrc, .zshrc, .cshrc).
#
# NOTE:
# - We use <PKG>_ROOT_TEMP variables because all variables are exported to 
#	environment when this script is source'd. In some cases, this may cause 
#	conflicts (e.g. on Pleiades, we use the module copy of PETSc, which 
#	defines PETSC_ROOT).

# Silence `zsh: no matches found: <file>`
if [[ -n "$ZSH_VERSION" ]]; then
	setopt +o nomatch 1> /dev/null 2>&1
fi

## Functions
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

#############################
# Compilers / runtime / SDKs
#############################
EMSCRIPTEN_ROOT_TEMP="${ISSM_EXT_DIR}/emscripten/install"
if [ -d ${EMSCRIPTEN_ROOT_TEMP} ]; then
	export EMSCRIPTEN_ROOT="${EMSCRIPTEN_ROOT_TEMP}" # Used in JavaScript build in installation of GSL, Triangle
fi

#############################
# Build systems
#############################
AUTOTOOLS_ROOT_TEMP="${ISSM_EXT_DIR}/autotools/install"
if [ -d "${AUTOTOOLS_ROOT_TEMP}" ]; then
	path_prepend "${AUTOTOOLS_ROOT_TEMP}/bin"
fi

CMAKE_ROOT_TEMP="${ISSM_EXT_DIR}/cmake/install"
if [ -d "${CMAKE_ROOT_TEMP}" ]; then
	path_prepend "${CMAKE_ROOT_TEMP}/bin"
fi

#############################
# Libraries / binaries
#############################

# NOTE: The following checks *must* come before PETSc as we prefer packages 
#		installed via PETSc
#
HDF5_ROOT_TEMP="${ISSM_EXT_DIR}/hdf5/install"
if [ -d "${HDF5_ROOT_TEMP}" ]; then
	export HDF5_ROOT="${HDF5_ROOT_TEMP}" # Used in installation of NetCDF, GDAL
	cpath_append "${HDF5_ROOT_TEMP}/include"
	library_path_prepend "${HDF5_ROOT_TEMP}/lib"
	dyld_library_path_prepend "${HDF5_ROOT_TEMP}/lib"
	ld_library_path_prepend "${HDF5_ROOT_TEMP}/lib"
fi

MUMPS_ROOT_TEMP="${ISSM_EXT_DIR}/mumps/install"
if [ -d "${MUMPS_ROOT_TEMP}" ]; then
	export MUMPS_ROOT="${MUMPS_ROOT_TEMP}" # Used in installation of PETSc
	library_path_append "${MUMPS_ROOT_TEMP}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${MUMPS_ROOT_TEMP}/lib" # Allows dynamic loader to find DLLs
	fi
fi

ZLIB_ROOT_TEMP="${ISSM_EXT_DIR}/zlib/install"
if [ -d "${ZLIB_ROOT_TEMP}" ]; then
	export ZLIB_ROOT="${ZLIB_ROOT_TEMP}" # Used in installation of NetCDF, GDAL, GMT
	ld_library_path_append "${ZLIB_ROOT_TEMP}/lib"
fi

PETSC_ROOT_TEMP="${ISSM_EXT_DIR}/petsc/install"
if [ -d "${PETSC_ROOT_TEMP}" ]; then
	export PETSC_ROOT="${PETSC_ROOT_TEMP}" # Used in installation of Gmsh
	cpath_prepend "${PETSC_ROOT_TEMP}/include"
	library_path_prepend "${PETSC_ROOT_TEMP}/lib"
	ld_library_path_prepend "${PETSC_ROOT_TEMP}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${PETSC_ROOT_TEMP}/lib" # Allows dynamic loader to find DLLs
	fi

	# In case we have installed certain external packages via PETSc
	#

	# BLAS
	if ls ${PETSC_ROOT_TEMP}/lib/libblas.* 1> /dev/null 2>&1 || ls ${PETSC_ROOT_TEMP}/lib/libfblas.* 1> /dev/null 2>&1; then
		export BLAS_ROOT="${PETSC_ROOT_TEMP}" # Used in installation of Dakota, GMT
	fi

	# HDF5
	if ls ${PETSC_ROOT_TEMP}/lib/libhdf5.* 1> /dev/null 2>&1; then
		export HDF5_ROOT="${PETSC_ROOT_TEMP}" # Used in installation of NetCDF, GDAL
		cpath_append "${PETSC_ROOT_TEMP}/include"
		library_path_append "${PETSC_ROOT_TEMP}/lib"
		dyld_library_path_append "${PETSC_ROOT_TEMP}/lib"
		ld_library_path_append "${PETSC_ROOT_TEMP}/lib"
	fi

	# LAPACK
	if ls ${PETSC_ROOT_TEMP}/lib/liblapack.* 1> /dev/null 2>&1 || ls ${PETSC_ROOT_TEMP}/lib/libflapack.* 1> /dev/null 2>&1; then
		export LAPACK_ROOT="${PETSC_ROOT_TEMP}" # Used in installation of Dakota, GMT
	fi

	# METIS
	if ls ${PETSC_ROOT_TEMP}/lib/libmetis.* 1> /dev/null 2>&1; then
		export METIS_ROOT="${PETSC_ROOT_TEMP}" # Used in installation of Gmsh
	fi

	# MPICH
	if [ -f "${PETSC_ROOT_TEMP}/bin/mpiexec" ]; then
		export MPI_ROOT=${PETSC_ROOT_TEMP}
		export MPI_DIR=${PETSC_ROOT_TEMP}
		export MPI_HOME=${PETSC_ROOT_TEMP} # Used in installation of Dakota
		export MPI_INC_DIR="${PETSC_ROOT_TEMP}/include"
		path_prepend "${PETSC_ROOT_TEMP}/bin"
		cpath_prepend "${PETSC_ROOT_TEMP}/include"
		ld_library_path_append "${PETSC_ROOT_TEMP}/lib"
	fi

	# ZLIB
	if ls ${PETSC_ROOT_TEMP}/lib/libz.* 1> /dev/null 2>&1; then
		export ZLIB_ROOT="${PETSC_ROOT_TEMP}" # Used in installation of NetCDF, GDAL
		ld_library_path_append "${PETSC_ROOT_TEMP}/lib"
	fi
fi

MPLAPACK_ROOT_TEMP="${ISSM_EXT_DIR}/mplapack/install"
if [ -d "${MPLAPACK_ROOT_TEMP}" ]; then
	cplus_include_path_prepend "${MPLAPACK_ROOT_TEMP}/include"
	cplus_include_path_prepend "${MPLAPACK_ROOT_TEMP}/include/mplapack"
	library_path_prepend "${MPLAPACK_ROOT_TEMP}/lib"
	ld_library_path_prepend "${MPLAPACK_ROOT_TEMP}/lib"
fi

ADJOINTPETSC_TEMP="${ISSM_EXT_DIR}/adjointpetsc/install"
if [ -d "${ADJOINTPETSC_TEMP}" ]; then
	export ADJOINTPETSC="${ADJOINTPETSC_TEMP}"
	ld_library_path_append "${ADJOINTPETSC_TEMP}/lib"
fi

SCOTCH_ROOT_TEMP="${ISSM_EXT_DIR}/scotch/install"
ld_library_path_append "${SCOTCH_ROOT_TEMP}/lib"

BOOST_ROOT_TEMP="${ISSM_EXT_DIR}/boost/install"
if [ -d "${BOOST_ROOT_TEMP}" ]; then
	export BOOST_ROOT=${BOOST_ROOT_TEMP} # Used in installation of Dakota
	export BOOST_DIR=${BOOST_ROOT_TEMP}
	export BOOSTROOT=${BOOST_ROOT_TEMP}
	path_append "${BOOST_ROOT_TEMP}/bin"
	library_path_prepend "${BOOST_ROOT_TEMP}/lib"
	ld_library_path_prepend "${BOOST_ROOT_TEMP}/lib"
	dyld_library_path_prepend "${BOOST_ROOT_TEMP}/lib"
fi

DAKOTA_ROOT_TEMP="${ISSM_EXT_DIR}/dakota/install"
if [ -d "${DAKOTA_ROOT_TEMP}" ]; then
	path_append "${DAKOTA_ROOT_TEMP}/bin"
	ld_library_path_prepend "${DAKOTA_ROOT_TEMP}/lib"
	dyld_library_path_prepend "${DAKOTA_ROOT_TEMP}/lib"
fi

CPPCHECK_ROOT_TEMP="${ISSM_EXT_DIR}/cppcheck/install"
path_append "${CPPCHECK_ROOT_TEMP}/bin"

GSL_ROOT_TEMP="${ISSM_EXT_DIR}/gsl/install"
if [ -d "${GSL_ROOT_TEMP}" ]; then
	export GSL_HOME="${GSL_ROOT_TEMP}" # Used in installation of Dakota
	cpath_prepend "${GSL_ROOT_TEMP}/include"
	ld_library_path_append "${GSL_ROOT_TEMP}/lib"
fi

NETCDF_ROOT_TEMP="${ISSM_EXT_DIR}/netcdf/install"
if [ -d "${NETCDF_ROOT_TEMP}" ]; then
	export NETCDF_ROOT="${NETCDF_ROOT_TEMP}" # Used in installation of GDAL, GMT
	path_prepend "${NETCDF_ROOT_TEMP}/bin"
	cpath_prepend "${NETCDF_ROOT_TEMP}/include"
	library_path_prepend "${NETCDF_ROOT_TEMP}/lib"
	ld_library_path_prepend "${NETCDF_ROOT_TEMP}/lib"
	dyld_library_path_prepend "${NETCDF_ROOT_TEMP}/lib"
fi

CURL_ROOT_TEMP="${ISSM_EXT_DIR}/curl/install"
if [ -d "${CURL_ROOT_TEMP}" ]; then
	export CURL_ROOT="${CURL_ROOT_TEMP}" # Used in installation of NetCDF, GDAL, GMT
	cpath_prepend "${CURL_ROOT_TEMP}/include"
	ld_library_path_prepend "${CURL_ROOT_TEMP}/lib"
	dyld_library_path_prepend "${CURL_ROOT_TEMP}/lib"
	path_prepend "${CURL_ROOT_TEMP}/bin"
fi

SQLITE_ROOT_TEMP="${ISSM_EXT_DIR}/sqlite/install"
if [ -d "${SQLITE_ROOT_TEMP}" ]; then
	export SQLITE_ROOT="${SQLITE_ROOT_TEMP}" # Used in installation of GDAL
	path_prepend "${SQLITE_ROOT_TEMP}/bin"
	cpath_prepend "${SQLITE_ROOT_TEMP}/include"
	library_path_prepend "${SQLITE_ROOT_TEMP}/lib"
	ld_library_path_prepend "${SQLITE_ROOT_TEMP}/lib"
fi

PROJ_ROOT_TEMP="${ISSM_EXT_DIR}/proj/install"
if [ -d "${PROJ_ROOT_TEMP}" ]; then
	export PROJ_ROOT="${PROJ_ROOT_TEMP}" # Used in installation of GDAL
	export PROJ_DATA="${PROJ_ROOT_TEMP}/share/proj" # In order to find proj.db (After PROJ 9.1)
	export PROJ_LIB="${PROJ_ROOT_TEMP}/share/proj"  # In order to find proj.db (Prior to PROJ 9.1)
	path_prepend "${PROJ_ROOT_TEMP}/bin"
	ld_library_path_prepend "${PROJ_ROOT_TEMP}/lib"
	dyld_library_path_prepend "${PROJ_ROOT_TEMP}/lib"
fi

GDAL_ROOT_TEMP="${ISSM_EXT_DIR}/gdal/install"
if [ -d "${GDAL_ROOT_TEMP}" ]; then
	export GDAL_ROOT="${GDAL_ROOT_TEMP}" # Used in installation of GMT
	path_prepend "${GDAL_ROOT_TEMP}/bin"
	library_path_prepend "${GDAL_ROOT_TEMP}/lib"
	ld_library_path_prepend "${GDAL_ROOT_TEMP}/lib"
	dyld_library_path_prepend "${GDAL_ROOT_TEMP}/lib"
fi

GSHHG_ROOT_TEMP="${ISSM_EXT_DIR}/gshhg/install"
if [ -d "${GSHHG_ROOT_TEMP}" ]; then
	export GSHHG_ROOT="${GSHHG_ROOT_TEMP}" # Used in installation of GMT
fi

GMT_ROOT_TEMP="${ISSM_EXT_DIR}/gmt/install"
if [ -d "${GMT_ROOT_TEMP}" ]; then
	path_prepend "${GMT_ROOT_TEMP}/bin"
	ld_library_path_append "${GMT_ROOT_TEMP}/lib"
	dyld_library_path_append "${GMT_ROOT_TEMP}/lib"
fi

GMSH_ROOT_TEMP="${ISSM_EXT_DIR}/gmsh/install"
if [ -d "${GMSH_ROOT_TEMP}" ]; then
	path_prepend "${GMSH_ROOT_TEMP}/bin"
	ld_library_path_append "${GMSH_ROOT_TEMP}/lib"
	dyld_library_path_append "${GMSH_ROOT_TEMP}/lib"
fi

TRIANGLE_ROOT_TEMP="${ISSM_EXT_DIR}/triangle/install"
if [ -d "${TRIANGLE_ROOT_TEMP}" ]; then
	ld_library_path_append "${TRIANGLE_ROOT_TEMP}/lib"
	dyld_library_path_append "${TRIANGLE_ROOT_TEMP}/lib"

	if [[ ${MINGW} -eq 1 ]]; then
		path_append "${TRIANGLE_ROOT_TEMP}/lib" # Allows dynamic loader to find DLLs
	fi
fi

BBFTP_ROOT_TEMP="${ISSM_EXT_DIR}/bbftp/install"
if [ -d "${BBFTP_ROOT_TEMP}" ]; then
	path_append "${BBFTP_ROOT_TEMP}/bin"
fi

SHAPELIB_ROOT_TEMP="${ISSM_EXT_DIR}/shapelib/install"
if [ -d "${SHAPELIB_ROOT_TEMP}" ]; then
	path_append "${SHAPELIB_ROOT_TEMP}/exec"
fi

ESMF_ROOT_TEMP="${ISSM_EXT_DIR}/esmf/install"
if [ -d "${ESMF_ROOT_TEMP}" ]; then
	path_prepend "${ESMF_ROOT_TEMP}/bin"
	ld_library_path_append "${ESMF_ROOT_TEMP}/lib/libO/Linux.gfortran.64.mpich.default"
fi

NEOPZ_ROOT_TEMP="${ISSM_EXT_DIR}/neopz/install"
if [ -d "${NEOPZ_ROOT_TEMP}" ]; then
	export REFPATTERNDIR="${NEOPZ_ROOT_TEMP}/include/refpatterns"
fi

YAMS_ROOT_TEMP="${ISSM_EXT_DIR}/yams/install"
if [ -d "${YAMS_ROOT_TEMP}" ]; then
	path_prepend "${YAMS_ROOT_TEMP}"
fi

VALGRIND_ROOT_TEMP="${ISSM_DIR}/externalpackages/valgrind/install"
if [ -d "${VALGRIND_ROOT_TEMP}" ]; then
	path_prepend "${VALGRIND_ROOT_TEMP}/bin"
fi

SHELL2JUNIT_ROOT_TEMP="${ISSM_EXT_DIR}/shell2junit/install"
if [ -d "${SHELL2JUNIT_ROOT_TEMP}" ]; then
	path_append "${SHELL2JUNIT_ROOT_TEMP}"
fi
