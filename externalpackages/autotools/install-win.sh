#!/bin/bash
set -eu


## Constants
#
AUTOCONF_VER="2.69"
AUTOMAKE_MIN_VER="1.16"
AUTOMAKE_VER="${AUTOMAKE_MIN_VER}.1"
LIBTOOL_VER="2.4.6"
M4_VER="1.4.19"

PREFIX="${ISSM_DIR}/externalpackages/autotools/install" # Set to location where external package should be installed

## Environment
#
export PATH="${PREFIX}/bin:${PATH}"

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX}

# Install m4
echo " === INSTALLING M4 ==="
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/m4-${M4_VER}.tar.gz" "m4-${M4_VER}.tar.gz"
tar -zxvf m4-${M4_VER}.tar.gz
mv m4-${M4_VER} src
cd src

./configure --prefix="${PREFIX}"
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
cd ..

# Install Autoconf
echo " === INSTALLING AUTOCONF =="
rm -rf src
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/autoconf-${AUTOCONF_VER}.tar.gz" "autoconf-${AUTOCONF_VER}.tar.gz"
tar -zxvf autoconf-${AUTOCONF_VER}.tar.gz
mv autoconf-${AUTOCONF_VER} src
cd src
./configure --prefix="${PREFIX}"
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
cd ..

# Install Automake
echo " === INSTALLING AUTOMAKE ==="
rm -rf src
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/automake-${AUTOMAKE_VER}.tar.gz" "automake-${AUTOMAKE_VER}.tar.gz"
tar -zxvf automake-${AUTOMAKE_VER}.tar.gz
mv automake-${AUTOMAKE_VER} src
cd src
./configure --prefix="${PREFIX}"
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
cd ..

# Install libtool
echo " === INSTALLING LIBTOOL ==="
rm -rf src
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/libtool-${LIBTOOL_VER}.tar.gz" "libtool-${LIBTOOL_VER}.tar.gz"
tar -zxvf libtool-${LIBTOOL_VER}.tar.gz
mv libtool-${LIBTOOL_VER} src
cd src
./configure --prefix="${PREFIX}"
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
cd ..

# This patch takes care of removing options passed to the linker that causes
# the build to fail, as well as changing some flags to match up to Microsoft
# compilers.
patch ./install/share/aclocal/libtool.m4 < ./patches/libtool.m4.patch

# Small change to Automake in order to get the right flags for Microsoft's
# compiler.
patch ./install/bin/automake < ./patches/automake.patch

# This patch is for ar-lib, and removes carriage return characters that cause
# commands to overwrite themselves and be misinterpreted during linking on
# Cygwin Windows.
patch ./install/share/automake-${AUTOMAKE_MIN_VER}/ar-lib < ./patches/ar-lib.patch
