##############################################################################
#
# Template CMake Configuration File.
#
##############################################################################
# The following CMake variables represent the minimum set of variables
# that are required to allow Dakota to
#   * find all prerequisite third party libraries (TPLs)
#   * configure compiler and MPI options
#   * set Dakota install path
#
# Instructions:
# 1. Read Dakota/INSTALL - Source Quick Start to use this template file.
#
# 2. Uncomment CMake variables below ONLY for values you need to change for
#    your platform. Edit variables as needed.
#
#    For example, if you are using a custom install of Boost, installed in
#    /home/me/usr/boost, uncomment both CMake Boost variables  and edit
#    paths:
#       set(BOOST_ROOT
#           "/home/me/usr/boost"
#           CACHE PATH "Use non-standard Boost install" FORCE)
#       set( Boost_NO_SYSTEM_PATHS TRUE
#            CACHE BOOL "Supress search paths other than BOOST_ROOT" FORCE)
#
#    Save file and exit.
#
# 6. Run CMake with script file. At terminal window, type:
#      $ cmake -C BuildCustom.cmake $DAK_SRC
#
#    If you have not followed instructions in INSTALL -Source Quick Start,
#    you will need to replace BuildCustom.cmake with the actual filename of
#    this file and $DAK_SRC with the actual path to Dakota source.
#
##############################################################################

##############################################################################
# Set BLAS, LAPACK library paths ONLY if in non-standard locations
##############################################################################
set( BLAS_LIBS
      "$ENV{BLAS_LIBS}"
      CACHE FILEPATH "Use non-standard BLAS library path" FORCE )
set( LAPACK_LIBS
      "$ENV{LAPACK_LIBS}"
      CACHE FILEPATH "Use non-standard BLAS library path" FORCE )

##############################################################################
# Set additional compiler options
# Uncomment and replace <flag> with actual compiler flag, e.g. -xxe4.2
##############################################################################
#set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} <flag>" 
#     CACHE STRING "C Flags my platform" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_IGNORE_CXX_SEEK"
     CACHE STRING "CXX Flags for my platform" )
#set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} <flag>"
#     CACHE STRING "Fortran Flags for my platform" )

##############################################################################
# Set MPI options
# Recommended practice is to set DAKOTA_HAVE_MPI and set MPI_CXX_COMPILER
# to a compiler wrapper.
##############################################################################
set( DAKOTA_HAVE_MPI ON
     CACHE BOOL "Build with MPI enabled" FORCE)
#set( MPI_CXX_COMPILER "path/to/mpicxx"
#     CACHE FILEPATH "Use MPI compiler wrapper" FORCE)

##############################################################################
# Set Boost path if CMake cannot find your installed version of Boost or
# if you have a custom Boost install location.
##############################################################################
set(BOOST_ROOT
    $ENV{BOOST_ROOT}
    CACHE PATH "Use non-standard Boost install" FORCE)
set( Boost_NO_SYSTEM_PATHS TRUE
     CACHE BOOL "Supress search paths other than BOOST_ROOT" FORCE)

##############################################################################
# Set Trilinos path if you want have a custom Trilinos install location. If
# not set, the Trilinos package, teuchos, will be build during the Dakota
# build.
##############################################################################
#set( Trilinos_DIR
#      "path/to/Trilinos/install"
#      CACHE PATH "Path to installed Trilinos" FORCE )

##############################################################################
# Customize DAKOTA
##############################################################################
set( CMAKE_INSTALL_PREFIX
     $ENV{DAK_INSTALL}
     CACHE PATH "Path to Dakota installation" )
