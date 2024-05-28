#!/bin/bash
set -eu

#Some cleanup
rm -rf install src
rm -rf vim72
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/vim-7.2.tar.bz2' 'vim-7.2.tar.bz2'

#Untar 
bzip2 -d -k vim-7.2.tar.bz2
tar -xvf vim-7.2.tar
rm vim-7.2.tar

#Move vim into install directory
mv vim72/* src
rm -rf vim72

#Configure vim (icc seems to have issues with wctype.h)
export CC=gcc
cd src/src 
./configure \
	--prefix="$ISSM_DIR/externalpackages/vim/install" \
	--with-gcc="/usr/bin/gcc" \
	--with-tlib="/lib/"

#Compile vim
make
make  install
