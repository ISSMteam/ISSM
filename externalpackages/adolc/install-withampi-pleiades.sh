#!/bin/bash
set -eu

#Some cleanup
rm -rf install src

# Keeping the following commented line for potential future use.
#git clone https://gitlab.com/adol-c/adol-c.git src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/ADOL-C.tar.gz' 'ADOL-C.tar.gz'

#Untar ADOL-C
tar -zxf  ADOL-C.tar.gz

#Compile ADOL-C
cd src
./configure --prefix=$ISSM_DIR/externalpackages/adolc/install  \
	--libdir=$ISSM_DIR/externalpackages/adolc/install/lib \
	--with-mpi-root="/nasa/sgi/mpt/2.06rp16/" \
	--enable-ampi \
	--with-ampi=$ISSM_DIR/externalpackages/adjoinablempi/install \
	--with-soname=adolc \
	--disable-tapedoc-values

make clean
if [ $# -eq 0 ]; then
	make V=1
else
	make -j $1 V=1
fi
make V=1 install
