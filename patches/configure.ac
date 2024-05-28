# Process this file with autoconf to produce a configure script.

#AUTOCONF
AC_INIT([ISSM],[4.2.5],[issm@jpl.nasa.gov],[issm],[http://issm.jpl.nasa.gov]) #Initializing configure
AC_CONFIG_AUX_DIR([./aux-config])         #Put config files in aux-config
AC_CONFIG_MACRO_DIR([m4])                 #m4 macros are located in m4
m4_include([m4/issm_options.m4])
AC_CANONICAL_TARGET                       #Determine the system type

#Compilers
AC_PROG_CC([icc cccl cl icl gcc])
AC_PROG_CPP
AC_PROG_CXX([icpc cccl cl icl g++])
AC_PROG_F77([ifort g77 gfortran])

#Libraries and linking
AC_F77_LIBRARY_LDFLAGS
AC_PATH_XTRA      #figure out X library and include paths
AC_PROG_RANLIB

#AUTOMAKE
#Initialize automake and declare foreign so that we don't need a ChangeLog, INSTALL, etc
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS([./config.h])   #Config file must be config.h
AM_SILENT_RULES([yes])           #Do not show compilation command by default
AM_PROG_CC_C_O
AM_PROG_AR

#Libtool
#LT_INIT

#Run issm_options.m4
ISSM_OPTIONS

#List all Makefiles
AC_CONFIG_FILES([Makefile
			src/Makefile
			src/c/Makefile 
			src/wrappers/Makefile
			src/wrappers/matlab/Makefile
			src/m/Makefile])

#End of configure.ac
AC_OUTPUT
