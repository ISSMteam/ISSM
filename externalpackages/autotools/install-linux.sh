#!/bin/bash
set -eu


## Constants
#
AUTOCONF_VER="2.69"
AUTOMAKE_VER="1.16.1"
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

## Fixes required by glibc-2.28
#
# Source: http://www.linuxfromscratch.org/lfs/view/development/chapter06/m4.html
#
sed -i 's/IO_ftrylockfile/IO_EOF_SEEN/' lib/*.c
echo "#define _IO_IN_BACKUP 0x100" >> lib/stdio-impl.h

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
echo " === INSTALLING AUTOCONF ==="
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
