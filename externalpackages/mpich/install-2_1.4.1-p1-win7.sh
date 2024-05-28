#!/bin/bash

#In order to install MPICH2 on your platform, go to www.mpich.org/downloads 
#and download the following file (Windows x86_64 for the Windows distribution): 
#mpich2-1.4.1p1-win-x86-64.msi. This file is also hosted on the ISSM website.


#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/mpich2-1.4.1p1-win-x86-64.msi' 'mpich2-1.4.1p1-win-x86-64.msi'

#once installed, create a symbolic link between the MPICH2 directory 
#and the install directory. For example: 
#ln -s /cygdrive/c/Program\ Files/MPICH2 install
