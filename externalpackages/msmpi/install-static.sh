#!/bin/bash
#
# Creates a local MS-MPI directory to be used for including MS-MPI with ISSM 
# distributable for Windows.
#
# Assumes that both MSYS2 pacman distribution of MS-MPI for the UCRT 
# environment (SDK) and Microsoft distribution of MS-MPI (runtime) have been 
# installed.
#
# To install the MS-MPI SDK:
# - Start a MSYS2 UCRT64 terminal instance.
# - Run,
#
#       pacman -S --needed mingw-w64-ucrt-x86_64-msmpi
#
# - Note the version of MS-MPI that was installed by running,
#
#       pacman -Qi mingw-w64-ucrt-x86_64-msmpi
#
#   For example, if this outputs '10.1.1-27', then the version is '10.1.1'.
#
# To install the MS-MPI runtime:
# - Navigate to the 'Microsoft Download Center' webpage: 
#   https://www.microsoft.com/en-us/download
# - Search for 'MPI <VERSION>' where <VERSION> is the version of the MS-MPI SDK 
#   installed previously via pacman.
# - Click the 'Download' button.
# - Only check the box for `msmpisetup.exe`, then click the 'Download' button.
# - When the download completes, run the installer.
# - Follow the prompts, using the default installation directory.
#


## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/msmpi/install"

MSMPI_RUNTIME_DIR=$(cygpath -u $(cygpath -ms "/c/Program Files/Microsoft MPI")
MSMPI_SDK_DIR="/ucrt64"

MSMPI_LIB="/c/Windows/System32/msmpi.dll"
MSMPIFEC_LIB="${MSMPI_SDK_DIR}/bin/libmsmpifec.dll"

# Cleanup
rm -rf ${PREFIX}
mkdir -p ${PREFIX} ${PREFIX}/bin ${PREFIX}/licenses

# # Generate static copy of MS-MPI library
# gendef ${MSMPI_LIB}
# dlltool -d msmpi.def -l libmsmpi.a -D ${MSMPI_LIB}

# # Copy MS-MPI library to 'lib' directory
# cp libmsmpi.a ${PREFIX}/lib
# cp libmsmpi.a ${PREFIX}/lib/msmpi.a
# cp libmsmpi.a ${PREFIX}/lib/msmpi.lib

# Copy MS-MPI library to distribution directory
cp ${MSMPI_LIB} ${PREFIX}/bin

# # Create link to shared library so that libtool can find it
# cd ${PREFIX}/lib
# ln -s msmpi.dll libmsmpi.dll

# Copy MS-MPI runtime files to respective distribution directories
cp -R ${MSMPI_RUNTIME_DIR}/Bin/* ${PREFIX}/bin
cp -R ${MSMPI_RUNTIME_DIR}/License/* ${PREFIX}/bin # For now, just store license files alongside binaries

# Copy MS-MPI SDK files to respective distribution directories
cp ${MSMPIFEC_LIB} ${PREFIX}/bin

