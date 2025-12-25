#!/bin/bash
set -eu

## Sources
# - https://earthsystemmodeling.org/docs/nightly/develop/ESMC_crefdoc/node5.html#SECTION05063000000000000000
# - https://cpp.hotexamples.com/examples/-/-/ESMC_MeshAddNodes/cpp-esmc_meshaddnodes-function-examples.html

## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/math77/install" # Set to location where external package should be installed

VER="8_0_1"

export ESMF_DIR="${ISSM_DIR}/externalpackages/esmf/src"
export ESMF_INSTALL_PREFIX="${PREFIX}"

# Cleanup
rm -rf ${PREFIX} src

# Download source
${ISSM_DIR}/scripts/DownloadExternalPackage.sh "https://github.com/esmf-org/esmf/archive/refs/tags/ESMF_${VER}.tar.gz" "ESMF_${VER}.tar.gz"

# Unpack source
tar -zxvf ESMF_${VER}.tar.gz
mv esmf-ESMF_${VER} ${ESMF_DIR}

# Compile and install
cd ${ESMF_DIR}
make
make install
