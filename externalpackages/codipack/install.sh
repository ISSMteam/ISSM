#!/bin/bash
set -eu

VER=1.8

# Cleanup from previous installation
rm -rf install CoDiPack-$VER.tar.gz

#Download development version
git clone https://github.com/SciCompKL/CoDiPack.git install

## Download source
#$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/CoDiPack-${VER}.tar.gz" "CoDiPack-${VER}.tar.gz"
#
## Untar
#tar -zxvf CoDiPack-$VER.tar.gz
#
## Move source into install directory
#mv CoDiPack-$VER install
#rm -rf CoDiPack-$VER/
