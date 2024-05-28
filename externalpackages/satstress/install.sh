#!/bin/bash
set -eu

#Some cleanup
rm -rf install SatStress-0.1.2
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/SatStress-0.1.2.tar.gz' 'SatStress-0.1.2.tar.gz'

#Untar 
tar -zxvf  SatStress-0.1.2.tar.gz

#Move SatStress into install directory
mv SatStress-0.1.2/* install
rm -rf SatStress-0.1.2

#Compile SatStress
cd install 
make test
