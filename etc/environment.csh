# Modifies path-related envrionment variables based on which external packages
# have been installed.
#
# ISSM_DIR and ISSM_ARCH should have been defined already in your shell
# settings file (i.e. .bashrc, .cshrc).
#
# TODO:
# - Condition all path modifications on existence of external package 'install'
#	directory
#

# Silence `<command>: no match`
set nonomatch=1

setenv PATH "{$PATH}:{$ISSM_DIR}/aux-config"
setenv PATH "{$PATH}:{$ISSM_DIR}/scripts"

set ISSM_EXT_DIR="{$ISSM_DIR}/externalpackages" # Redefine this constant if externalpackages are installed to a different directory

#############################
# Build systems
#############################
set AUTOTOOLS_ROOT="{$ISSM_EXT_DIR}/autotools/install"
setenv PATH "{$AUTOTOOLS_ROOT}/bin:{$PATH}"

set CMAKE_ROOT="{$ISSM_EXT_DIR}/cmake/install"
setenv PATH "{$CMAKE_ROOT}/bin:{$PATH}"

#############################
# Libraries / binaries
#############################
set MPI_ROOT="{$ISSM_EXT_DIR}/mpich/install"
if ( -d {$MPI_ROOT} ) then
	setenv MPI_DIR {$MPI_ROOT}
	setenv MPI_HOME {$MPI_ROOT} # Used in installation of Dakota
	setenv MPI_INC_DIR {$MPI_ROOT}/include
	setenv PATH "{$MPI_ROOT}/bin:{$PATH}"
	setenv CPATH "{$MPI_ROOT}/include:{$CPATH}"
	setenv LD_LIBRARY_PATH "{$LD_LIBRARY_PATH}:{$MPI_ROOT}/lib"
endif

# NOTE: Check *must* come before PETSc as we prefer packages installed via 
# 		PETSc
#
set ZLIB_ROOT="{$ISSM_EXT_DIR}/zlib/install"
if ( -d {$ZLIB_ROOT} ) then
	setenv ZLIB_ROOT {$ZLIB_ROOT} # Used in installation of NetCDF, GDAL, GMT
	setenv LD_LIBRARY_PATH "{$LD_LIBRARY_PATH}:{$LD_LIBRARY_PATH}/lib"
endif

set PETSC_ROOT="{$ISSM_EXT_DIR}/petsc/install"
if ( -d {$PETSC_ROOT} ) then
	setenv PETSC_ROOT {$PETSC_ROOT}
	setenv LD_LIBRARY_PATH {$PETSC_ROOT}/lib:{$LD_LIBRARY_PATH}

	# In case we have installed certain external packages via PETSc
	#

	# BLAS
	if ( `find {$PETSC_ROOT}/lib -name libblas.*` != "" || `find {$PETSC_ROOT}/lib -name libfblas.*` != "" ) then
		setenv BLAS_ROOT "{$PETSC_ROOT}" # Used in installation of Dakota, GMT
	endif

	# HDF5
	if ( `find {$PETSC_ROOT}/lib -name libhdf5.*` != "" ) then
		setenv HDF5_ROOT "{$PETSC_ROOT}" # Used in installation of NetCDF, GDAL
		setenv CPATH {$CPATH}:{$PETSC_ROOT}/include
		setenv LIBRARY_PATH {$LIBRARY_PATH}:{$PETSC_ROOT}/lib
		setenv DYLD_LIBRARY_PATH {$DYLD_LIBRARY_PATH}:{$PETSC_ROOT}/lib
		setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$PETSC_ROOT}/lib
	endif

	# LAPACK
	if ( `find {$PETSC_ROOT}/lib -name liblapack.*` != "" || `find {$PETSC_ROOT}/lib -name libflapack.*` != "" ) then
		setenv LAPACK_ROOT "{$PETSC_ROOT}" # Used in installation of Dakota, GMT
	endif

	# METIS
	if ( `find {$PETSC_ROOT}/lib -name libmetis.*` != "" ) then
		setenv METIS_ROOT "{$PETSC_ROOT}" # Used in installation of Gmsh
	endif

	# MPICH
	if ( -f "{$PETSC_ROOT}/bin/mpiexec" ) then
		set MPI_ROOT={$PETSC_ROOT}
		setenv MPI_DIR {$MPI_ROOT}
		setenv MPI_HOME {$MPI_ROOT} # Used in installation of Dakota
		setenv MPI_INC_DIR {$MPI_ROOT}/include
		setenv PATH {$MPI_ROOT}/bin:{$PATH}
		setenv CPATH {$MPI_ROOT}/include:{$CPATH}
		setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$MPI_ROOT}/lib
	endif

	# ZLIB
	if ( `find {$PETSC_ROOT}/lib -name libz.*` != "" ) then
		setenv ZLIB_ROOT "{$PETSC_ROOT}" # Used in installation of NetCDF, GDAL
		setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$PETSC_ROOT}/lib
	endif
endif

set DAKOTA_ROOT="{$ISSM_EXT_DIR}/dakota/install"
if ( -d {$DAKOTA_ROOT} ) then
	setenv PATH {$PATH}:{$DAKOTA_ROOT}/bin
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$DAKOTA_ROOT}/lib
	setenv DYLD_LIBRARY_PATH {$DAKOTA_ROOT}/lib:{$DYLD_LIBRARY_PATH}
endif

set BOOST_ROOT="{$ISSM_EXT_DIR}/boost/install"
if ( -d {$BOOST_ROOT} ) then
	setenv BOOST_ROOT {$BOOST_ROOT} # Used in installation of Dakota
	setenv BOOST_DIR {$BOOST_ROOT}
	setenv BOOSTROOT {$BOOST_ROOT}
	setenv LIBRARY_PATH {$BOOST_ROOT}/lib:{$LIBRARY_PATH}
	setenv LD_LIBRARY_PATH {$BOOST_ROOT}/lib:{$LD_LIBRARY_PATH}
	setenv DYLD_LIBRARY_PATH {$BOOST_ROOT}/lib:{$DYLD_LIBRARY_PATH}
	setenv PATH {$BOOST_ROOT}/bin:{$PATH}
endif

set GSL_ROOT="{$ISSM_EXT_DIR}/gsl/install"
if ( -d {$GSL_ROOT} ) then
	setenv GSL_HOME {$GSL_ROOT} # Used in installation of Dakota
	setenv GSL_ROOT {$GSL_ROOT}
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$GSL_ROOT}/lib
	setenv DYLD_LIBRARY_PATH {$GSL_ROOT}/lib:{$DYLD_LIBRARY_PATH}
endif

set NETCDF_ROOT="{$ISSM_EXT_DIR}/netcdf/install"
if ( -d {$NETCDF_ROOT} ) then
	setenv NETCDF_ROOT "{$NETCDF_ROOT}" # Used in installation of GDAL, GMT
	setenv PATH {$PATH}:{$NETCDF_ROOT}/bin
	setenv CPATH {$CPATH}:{$NETCDF_ROOT}/include
	setenv LIBRARY_PATH {$LIBRARY_PATH}:{$NETCDF_ROOT}/lib
	setenv DYLD_LIBRARY_PATH {$DYLD_LIBRARY_PATH}:{$NETCDF_ROOT}/lib
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$NETCDF_ROOT}/lib
endif

set CURL_ROOT="{$ISSM_EXT_DIR}/curl/install"
if ( -d {$CURL_ROOT} ) then
	setenv CURL_ROOT "{$CURL_ROOT}" # Used in installation of NetCDF, GDAL, GMT
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$CURL_ROOT}/lib
	setenv DYLD_LIBRARY_PATH {$DYLD_LIBRARY_PATH}:{$CURL_ROOT}/lib
	setenv PATH {$PATH}:{$CURL_ROOT}/bin
endif

set HDF5_ROOT="{$ISSM_EXT_DIR}/hdf5/install"
if ( -d {$HDF5_ROOT} ) then
	setenv HDF5_ROOT "{$HDF5_ROOT}" # Used in installation of NetCDF, GDAL
	setenv CPATH {$CPATH}:{$HDF5_ROOT}/include
	setenv LIBRARY_PATH {$LIBRARY_PATH}:{$HDF5_ROOT}/lib
	setenv DYLD_LIBRARY_PATH {$DYLD_LIBRARY_PATH}:{$HDF5_ROOT}/lib
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$HDF5_ROOT}/lib
endif

set SQLITE_ROOT="{$ISSM_EXT_DIR}/sqlite/install"
if ( -d {$SQLITE_ROOT} ) then
	setenv PATH {$PATH}:{$SQLITE_ROOT}/bin
	setenv LIBRARY_PATH {$LIBRARY_PATH}:{$SQLITE_ROOT}/lib
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$SQLITE_ROOT}/lib
endif

set PROJ_ROOT="{$ISSM_EXT_DIR}/proj/install"
if ( -d {$PROJ_ROOT} ) then
	setenv PROJ_ROOT "${PROJ_ROOT}" # Used in installation of GDAL
	setenv DYLD_LIBRARY_PATH {$DYLD_LIBRARY_PATH}:{$PROJ_ROOT}/lib
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$PROJ_ROOT}/lib
endif

set GDAL_ROOT="{$ISSM_EXT_DIR}/gdal/install"
if ( -d {$GDAL_ROOT} ) then
	setenv GDAL_ROOT "{$GDAL_ROOT}" # Used in installation of GMT
	setenv PATH {$GDAL_ROOT}/bin:{$PATH}
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$GDAL_ROOT}/lib
endif

set GSHHG_ROOT="{$ISSM_EXT_DIR}/gshhg/install"
if ( -d {$GSHHG_ROOT} ) then
	setenv GSHHG_ROOT "{$GSHHG_ROOT}" # Used in installation of GMT
endif

set GMT_ROOT="{$ISSM_EXT_DIR}/gmt/install"
if ( -d {$GMT_ROOT} ) then
	setenv PATH {$GMT_ROOT}/bin:{$PATH}
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$GMT_ROOT}/lib
	setenv DYLD_LIBRARY_PATH {$DYLD_LIBRARY_PATH}:{$GMT_ROOT}/lib
endif

set GMSH_ROOT="{$ISSM_EXT_DIR}/gmsh/install"
if ( -d {$GMSH_ROOT} ) then
	setenv PATH {$GMSH_ROOT}/bin:{$PATH}
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$GMSH_ROOT}/lib
	setenv DYLD_LIBRARY_PATH {$DYLD_LIBRARY_PATH}:{$GMSH_ROOT}/lib
endif

set TRIANGLE_ROOT="{$ISSM_EXT_DIR}/triangle/install"
if ( -d {$TRIANGLE_ROOT} ) then
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$TRIANGLE_ROOT}/lib
	setenv DYLD_LIBRARY_PATH {$DYLD_LIBRARY_PATH}:{$TRIANGLE_ROOT}/lib
endif

set VALGRIND_ROOT="{$ISSM_EXT_DIR}/valgrind/install"
if ( -d {$VALGRIND_ROOT} ) then
	setenv PATH {$VALGRIND_ROOT}/bin:{$PATH}
endif

set SHELL2JUNIT_ROOT="{$ISSM_EXT_DIR}/shell2junit/install"
if ( -d {$SHELL2JUNIT_ROOT} ) then
	setenv PATH {$PATH}:{$SHELL2JUNIT_ROOT}/install
endif
