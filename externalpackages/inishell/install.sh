#!/bin/bash
#inishell install package:  this package is not distributed by ISSM. Please request access to the code by 
#contacting Mathias Bavay () or logging onto http://models.slf.ch/ and registering.  Once you have a tarball 
#of the code, please use this script to install.

#we assume you have a inishell-src-*.tgz  tarball
set -eu

#Do some inquiry about the names of the tar balls: 
source_tar=`ls inishell-src-*.tgz`
source_version=`echo $source_tar | sed 's/inishell-src-//g' | sed 's/.tgz//g'`

if [[ $source_tar == "" ]]; then 
	echo "Could not identify a tar ball for the source code, name should be inishell-src-*.tgz"
	exit 1
fi
if [[ $source_version == "" ]]; then 
	echo "Could not identify a tar ball version for the source code"
	exit 1
fi


#Some cleanup
rm -rf src install inishell-$source_version
mkdir install

#First deal with source code 
tar -zxvf  $source_tar
mv inishell-$source_version src
rm -rf inishell-$source_version

#Build inishell
cd src && ant snowpack

#Put the .jar in the install directory
cp dist/inishell.jar ../install

#Install script to launch inishell jar file directly
cd ..
cp scripts/inishell install
