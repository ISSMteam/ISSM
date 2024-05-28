#!/bin/bash
#SNOWPACK install package:  this package is not distributed by ISSM. Please request access to the code by 
#contacting Mathias Bavay () or logging onto http://models.slf.ch/ and registering.  Once you have a tarball 
#of the code, please use this script to install.

#we assume you have a Snowpack-*-src.tar.gz  tarball
set -eu

#Do some inquiry about the names of the tar balls: 
source_tar=`ls Snowpack-*-src.tar.gz`
source_version=`echo $source_tar | sed 's/Snowpack-//g' | sed 's/-src.tar.gz//g'`

if [[ $source_tar == "" ]]; then 
	echo "Could not identify a tar ball for the source code, name should be Snowpack-*-src.tar.gz"
	exit 1
fi
if [[ $source_version == "" ]]; then 
	echo "Could not identify a tar ball version for the source code"
	exit 1
fi

#Some cleanup
rm -rf src install Snowpack-$source_version

#First deal with source code 
tar -zxvf  $source_tar
mv Snowpack-$source_version/usr src
rm -rf Snowpack-$source_version

#Build library
cd src && rm -rf CMakeCache.txt 
cmake -DCMAKE_INSTALL_PREFIX:PATH=$ISSM_DIR/externalpackages/snowpack/install \
	-DMETEOIO_INCLUDE_DIR:PATH=$ISSM_DIR/externalpackages/meteoio/install/include \
	-DMETEOIO_LIBRARY:PATH=$ISSM_DIR/externalpackages/meteoio/install/lib/libmeteoio.dylib\
	.
	#-DCMAKE_VERBOSE_MAKEFILE=true \

#Compile
if [ $# -eq 0 ]; then
	make  all install
else
	make -j $1  all install
fi

#Build binary: 
cd applications/snowpack
cmake -DCMAKE_MODULE_PATH:PATH=/Users/larour/issm-uci/trunk-jpl/externalpackages/snowpack/src/tools/cmake/      -DMETEOIO_INCLUDE_DIR:PATH=$ISSM_DIR/externalpackages/meteoio/install/include  .

#Compile
if [ $# -eq 0 ]; then
	make  all install
else
	make -j $1  all install
fi


