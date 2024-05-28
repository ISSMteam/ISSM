#!/bin/bash
set -eu 

#Some cleanup
rm -rf src install BLAS blas.tgz
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/blas.tgz' 'blas.tgz'

#Untar 
tar -zxvf  blas.tgz

#Move blas into install directory
mv BLAS/* src
rm -rf BLAS

#install
cd src 
cp ../configs/linux64/make.inc ./
make 

#Compile 
cd ../install
mkdir lib
cd lib
cp ../../src/*.a .
ln -s blas_LINUX.a blas.a
ln -s blas_LINUX.a libblas.a
