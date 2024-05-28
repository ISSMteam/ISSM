#!/bin/bash
set -eu

#Some cleanup
rm -rf install
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/gslib90_ls.tar.gz' 'gslib90_ls.tar.gz'

#Untar 
tar -zxvf  gslib90_ls.tar.gz

#Move gslib into install directory
mv gslib90/* install
rm -rf gslib90

#Change compiler to gfortran
cd install
cat Makefile | sed -e "s/FC=g95/FC=gfortran/g" > Makefile.bak
mv Makefile.bak Makefile
cat gslib/Makefile | sed -e "s/FC=g95/FC=gfortran/g" > Makefile.bak
mv Makefile.bak gslib/Makefile
make 
