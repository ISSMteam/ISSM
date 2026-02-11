
#Initializing configure
AC_INIT([Ice-sheet and Sea-level System Model (ISSM)],[2026.1],[https://github.com/ISSMteam/ISSM/],[issm],[https://issmteam.github.io/ISSM-Documentation/])
AC_CONFIG_AUX_DIR([./aux-config]) # Put config files in aux-config
AC_CONFIG_MACRO_DIR([m4])         # m4 macros are located in m4
m4_include([m4/issm_options.m4])

#print header
AC_MSG_NOTICE(============================================================================)
AC_MSG_NOTICE(=      AC_PACKAGE_NAME AC_PACKAGE_VERSION configuration      =)
AC_MSG_NOTICE(============================================================================)

#Determine System type and OS
AC_CANONICAL_TARGET

#Compilers
AC_PROG_CC( [mpicc icc cl icl gcc])             # C compiler
AC_PROG_CPP                                     # C preprocessor (after CC)
AC_PROG_CXX([mpicxx mpic++ icpc cl icl g++])    # C++ compiler
AC_PROG_F77([mpifort ifort g77 gfortran])       # Fortran 77 compiler
AC_PROG_FC( [mpifort ifort gfortran])           # Fortran compiler

#AUTOMAKE
AM_INIT_AUTOMAKE([foreign])     # Initialize automake and declare foreign (don't need a ChangeLog, INSTALL, etc)
AC_CONFIG_HEADERS([./config.h]) # Config file must be config.h
AM_SILENT_RULES([yes])          # Reduce compilation and linking output
AM_PROG_AR

#Libtool
LT_INIT([win32-dll])

# ---------------------------------------------------------------------
#  Run ISSM’s usual option-detection macro collection
# ---------------------------------------------------------------------
ISSM_OPTIONS

# ---------------------------------------------------------------------
#  Output files
# ---------------------------------------------------------------------
AC_CONFIG_FILES([
        Makefile
        src/Makefile
        src/c/Makefile
        src/m/Makefile
        src/wrappers/Makefile
        src/wrappers/python/Makefile
        src/wrappers/matlab/Makefile
        src/wrappers/javascript/Makefile
])

AC_OUTPUT
