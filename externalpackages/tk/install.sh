#!/bin/bash
set -eu

sudochoice=0;

#Some cleanup
rm -rf src
rm -rf install
rm -rf tk8.5.12
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/tk8.5.12-src.tar.gz' 'tk8.5.12.tar.gz'

#Untar 
tar -zxvf  tk8.5.12.tar.gz

#Move tk into src directory
mv tk8.5.12/* src
rm -rf tk8.5.12

cd src/unix

#User mode: 
if [[ $sudochoice == "0" ]]; 
then 
	./configure --prefix=$ISSM_DIR/externalpackages/tk/install
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
