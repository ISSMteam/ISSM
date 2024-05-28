#!/bin/bash
#METEOIO install package:  this package is not distributed by ISSM. Please request access to the code by 
#contacting Mathias Bavay () or logging onto http://models.slf.ch/ and registering.  Once you have a tarball 
#of the code, please use this script to install.

#we assume you have a MeteoIO-*-tar.gz  tarball
set -eu

#Do some inquiry about the names of the tar balls: 
source_tar=`ls MeteoIO-*.tar.gz`
source_version=`echo $source_tar | sed 's/MeteoIO-//g' | sed 's/.tar.gz//g'`

if [[ $source_tar == "" ]]; then 
	echo "Could not identify a tar ball for the source code, name should be MeteoIO-*-tar.gz"
	exit 1
fi
if [[ $source_version == "" ]]; then 
	echo "Could not identify a tar ball version for the source code"
	exit 1
fi

#Some cleanup
rm -rf src install MeteoIO-$source_version

#First deal with source code 
tar -zxvf  $source_tar
mv MeteoIO-$source_version/usr src
rm -rf MeteoIO-$source_version

#Reset makefile: 
cd src && rm -rf CMakeCache.txt 
cmake -DCMAKE_INSTALL_PREFIX:PATH=$ISSM_DIR/externalpackages/meteoio/install .

#Compile 
if [ $# -eq 0 ]; then
	make all install
else
	make -j $1 all install
fi
