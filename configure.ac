
#Initializing configure
AC_INIT([Ice-sheet and Sea-level System Model (ISSM)],[4.24],[https://issm.ess.uci.edu/forum/],[issm],[http://issm.jpl.nasa.gov])
AC_CONFIG_AUX_DIR([./aux-config])	# Put config files in aux-config
AC_CONFIG_MACRO_DIR([m4])			# m4 macros are located in m4
m4_include([m4/issm_options.m4])

#print header
AC_MSG_NOTICE(============================================================================)
AC_MSG_NOTICE(=      AC_PACKAGE_NAME AC_PACKAGE_VERSION configuration      =)
AC_MSG_NOTICE(============================================================================)

#Determine System type and OS
AC_CANONICAL_TARGET

#Compilers
AC_PROG_CC([icc cl icl gcc])
AC_PROG_CPP
AC_PROG_CXX([icpc cl icl g++])
AC_PROG_F77([ifort g77 gfortran])
AC_PROG_FC([ifort gfortran])

#AUTOMAKE
AM_INIT_AUTOMAKE([foreign])     # Initialize automake and declare foreign (don't need a ChangeLog, INSTALL, etc)
AC_CONFIG_HEADERS([./config.h]) # Config file must be config.h
AM_SILENT_RULES([yes])          # Reduce compilation and linking output
AM_PROG_AR

#Libtool
LT_INIT([win32-dll])

#Run issm_options.m4
ISSM_OPTIONS

#List all Makefiles
AC_CONFIG_FILES([
	Makefile
	src/Makefile
	src/c/Makefile
	src/wrappers/Makefile
	src/wrappers/python/Makefile
	src/wrappers/matlab/Makefile
	src/wrappers/javascript/Makefile
	src/m/Makefile
])

#End of configure.ac
AC_OUTPUT
