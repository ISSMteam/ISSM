#!/bin/bash
set -eu

#Some cleanup
rm -rf src install mpich-3.0.4
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/mpich-3.0.4.tar.gz' 'mpich-3.0.4.tar.gz'

#Untar 
tar -zxvf  mpich-3.0.4.tar.gz

#Move mpich into src directory
mv mpich-3.0.4/* src
rm -rf mpich-3.0.4

#Configure mpich
cd src
./configure \
	--prefix="$ISSM_DIR/externalpackages/mpich/install" \
	--enable-shared \
	CC=/usr/bin/gcc CXX=/usr/bin/g++ \
	F77=/usr/local/gfortran/bin/x86_64-apple-darwin10-gfortran \
	LDFLAGS="-L/usr/lib/ -lstdc++ -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin10/4.6.2/ -lgfortran"

	#CC=llvm-gcc \

#Compile mpich (this new version supports parallel make)
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
