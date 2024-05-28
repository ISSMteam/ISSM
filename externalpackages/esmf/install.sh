#!/bin/bash
set -eu

#DOC: https://earthsystemmodeling.org/docs/nightly/develop/ESMC_crefdoc/node5.html#SECTION05063000000000000000
# https://cpp.hotexamples.com/examples/-/-/ESMC_MeshAddNodes/cpp-esmc_meshaddnodes-function-examples.html

$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/ESMF_8_0_1.tar.gz" "ESMF_8_0_1.tar.gz"
tar -zxvf ESMF_8_0_1.tar.gz
mv ESMF_8_0_1 esmf
export ESMF_DIR=$ISSM_DIR/externalpackages/esmf/esmf
export ESMF_INSTALL_PREFIX=$ISSM_DIR/externalpackages/esmf/install

#Compile and install esmf
cd esmf
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
