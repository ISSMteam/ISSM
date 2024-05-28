#!/bin/bash
#
# Creates a local MS-MPI directory to be used for supplying MPI headers files 
# and libraries to ISSM configuration and certain external packages.
#
# Assumes that Microsoft MPI and MPI SDK have been installed. To do so,
# - Navigate to https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi
# - Under the 'MS-MPI Downloads' heading, click the link for 
#	'MS-MPI v<version>', where <version> is the latest version available (as of 
#	this writing, 10.1.2)
# - Click the 'Download' button
# - Make sure both boxes are checked
# - Click the 'Save File' button in each prompt
# - When the downloads are complete, run each installer
#
# TODO:
# - Commit Microsoft MPI and Microsoft SDK installers or source code to 
#	external packages source repository, then update this documentation to note 
#	that they are available
# - Attempt to download Microsoft MPI and Microsoft SDK installers or source 
#	code and (compile and) install with this script
# - Alternatively, instruct users to install MSYS2 MinGW 64-bit MS-MPI package 
#	with,
#
#		pacman -S mingw-w64-x86_64-msmpi
#
# remove this script, its parent directory, and references to it from 
# configuration files in $ISSM_DIR/jenkins directory and documentation
#


## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/msmpi/install"

MSMPI_BIN_DIR=$(cygpath -u $(cygpath -ms "/c/Program Files/Microsoft MPI/Bin"))
MSMPI_INC_DIR=$(cygpath -u $(cygpath -ms "/c/Program Files (x86)/Microsoft SDKs/MPI/Include"))
MSMPI_LIB="/c/Windows/System32/msmpi.dll"

# Cleanup
rm -rf ${PREFIX}
mkdir -p ${PREFIX} ${PREFIX}/bin ${PREFIX}/include ${PREFIX}/lib

# Copy MS-MPI header files to 'include' directory
cp ${MSMPI_INC_DIR}/mpi.h ${PREFIX}/include
cp ${MSMPI_INC_DIR}/mpi.f90 ${PREFIX}/include
cp ${MSMPI_INC_DIR}/mpif.h ${PREFIX}/include
cp ${MSMPI_INC_DIR}/mpio.h ${PREFIX}/include
cp ${MSMPI_INC_DIR}/x64/mpifptr.h ${PREFIX}/include

# Copy MS-MPI library to 'lib' directory
cp ${MSMPI_LIB} ${PREFIX}/lib

# Create link to shared library so that libtool can find it
cd ${PREFIX}/lib
ln -s msmpi.dll libmsmpi.dll
