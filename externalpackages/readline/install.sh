#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf readline-6.2.2
mkdir src 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/readline-6.2.2.tar.gz' 'readline-6.2.2.tar.gz'

#Untar 
tar -zxvf  readline-6.2.2.tar.gz

#Move readline into src directory
mv readline-6.2.2/* src
rm -rf readline-6.2.2

#install
cd src
python setup.py install
