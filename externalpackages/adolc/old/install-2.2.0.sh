#!/bin/bash
set -eu

#some issues on macosx64 with ISSM's autoconf. you might want to run native to mac on this.

#Some cleanup
rm -rf install ADOL-C-2.2.0 src trunk

#Create install directories
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/ADOL-C-2.2.0.tar.gz' 'ADOL-C-2.2.0.tar.gz'

#Untar 
tar -zxvf  ADOL-C-2.2.0.tar.gz

#Move ADOL-C into install directory
mv ADOL-C-2.2.0/* src
rm -rf ADOL-C-2.2.0

#Compile ADOL-C
cd src 

#export CC=gcc
#export CXX=g++
#export CFLAGS="-arch x86_64"
#export CXXFLAGS="-arch x86_64"

./configure \
	--prefix=$ISSM_DIR/externalpackages/adolc/install \
	--enable-sparse \
	--enable-docexa \
	--enable-addexa \
	--disable-shave

if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install


#Ok, bug with libtool: replace all LIBTOOL= by LIBTOOL=libtool 
#in all Makefiles
for i in `find ./ -name Makefile `
do
	echo $i
	cat $i | sed 's/LIBTOOL =/LIBTOOL = libtool/g' > $i.bak 
	mv $i.bak $i
done

#remake: 
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
