####################################################################
#  BLAS make include file.                                         #
#  March 2007                                                      #
####################################################################
#
SHELL = /usr/bin/bash
#
#  The machine (platform) identifier to append to the library names
#
PLAT =
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
FORTRAN  = gfortran
OPTS     = -O2 -fPIC
DRVOPTS  = $(OPTS)
NOOPT    = -O0 -fPIC
LOADER   = gfortran
LOADOPTS =
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib

MODULE			= blas
IMP_LIB_NAME	= lib$(MODULE).dll.a
LIB_NAME		= msys-$(MODULE).dll
#
#  The location and name of the Reference BLAS library.
#
BLASLIB = $(LIB_NAME)
