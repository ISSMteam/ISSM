#!/bin/bash
set -eu

sudochoice=0;

#Some cleanup
rm -rf src
rm -rf install
rm -rf tcl8.5.11
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/tcl8.5.11.tar.gz' 'tcl8.5.11.tar.gz'

#Untar 
tar -zxvf  tcl8.5.11.tar.gz

#Move tcl into src directory
mv tcl8.5.11/* src
rm -rf tcl8.5.11

#Configure tcl
ver="8.4.12"

cd src/unix

#User mode: 
if [[ $sudochoice == "0" ]]; 
then 
	./configure --prefix=$ISSM_DIR/externalpackages/tcl/install
	if [ $# -eq 0 ]; then
		make
	else
		make -j $1
	fi
	make install 
fi

#sudo mode: 
if [[ $sudochoice == "1" ]]; 
then
	sudo ./configure 
	if [ $# -eq 0 ]; then
		sudo make
	else
		sudo make -j $1
	fi
	sudo make install 
fi
