#!/bin/bash
set -eu

#clean up
rm -rf Cython-0.26

#download numpy first
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/Cython-0.26.tar.gz' 'cython.tar.gz'

#install numpy
tar -zxvf cython.tar.gz
cd Cython-0.26
python setup.py install
