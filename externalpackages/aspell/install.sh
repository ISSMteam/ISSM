#!/bin/bash
set -eu

#0: cleanup
#1:  install aspell
#2:  install en dictionary
step=2


#Some cleanup
if  [[ $step  == "0" ]]; then
	rm -rf src install dicts
	rm -rf aspell-0.50.5
	rm -rf aspell5-en-6.0.0
fi

#install aspell
if  [[ $step  == "1" ]]; then

	mkdir src install dicts

	#Download from ISSM server
	$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/aspell-0.50.5.tar.gz' 'aspell-0.50.5.tar.gz'

	#Untar 
	tar -zxvf  aspell-0.50.5.tar.gz

	#Move aspell into src directory
	mv aspell-0.50.5/* src
	rm -rf aspell-0.50.5

	#Configure aspell
	cd src
	./configure \
		--prefix="$ISSM_DIR/externalpackages/aspell/install" \

	#Compile and install aspell
	if [ -z $1 ]; then
		make
	else
		make -j $1
	fi
	make install
fi

#languages
if  [[ $step  == "2" ]]; then

	#Download from ISSM server
	$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/aspell5-en-6.0.0.tar.gz' 'aspell5-en-6.0.0.tar.gz'

	#Untar 
	tar -zxvf  aspell5-en-6.0.0.tar.gz

	#Move aspell into src directory
	mv aspell5-en-6.0.0 dicts

	#install
	cd dicts/aspell5-en-6.0.0
	./configure
	make install
fi
