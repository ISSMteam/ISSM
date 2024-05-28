#!/bin/bash

#Some cleanup
rm -rf install
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/yams2-linux.gz' 'yams2-linux.gz'
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/yams2-osx.gz' 'yams2-osx.gz'
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/yams2-win.exe' 'yams2-win.exe'

#loop over the binaries
for i in yams*
do
	name=$i;
	cp $i install/

	#uncompress if necessary
	if echo $i | grep -q ".gz"
	then
		gunzip install/$i
	fi

	#permissions
	chmod 777 install/*
done
