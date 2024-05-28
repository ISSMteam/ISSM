#!/bin/bash
set -eu
#modelE  downloaded from the gis repository of the GISS (Goddard Institute for Space Studies)
#at  http://www.giss.nasa.gov/tools/modelE/

#Some cleanup
rm -rf src install  modelE_AR5_branch
mkdir src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/modelE_AR5_branch.2012.03.13_10.12.21.tgz' 'modelE_AR5_branch.2012.03.13_10.12.21.tgz'

#Untar 
tar -zxvf  modelE_AR5_branch.2012.03.13_10.12.21.tgz

#Move modelE into install directory
mv modelE_AR5_branch/* src
rm -rf modelE_AR5_branch

#Apply patches
cd src 
