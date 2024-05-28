#!/bin/bash
set -eu

#Some cleanup
rm -rf install mercurial-1.7.3
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/mercurial-1.7.3.tar.gz' 'mercurial-1.7.3.tar.gz'

#Untar 
tar -zxvf  mercurial-1.7.3.tar.gz

#Move mercurial into install directory
mv mercurial-1.7.3/* install
rm -rf mercurial-1.7.3

#Apply patches
cd install 
#patch Lib/Makefile ../lib_Makefile.patch

#Compile mercurial
make local
