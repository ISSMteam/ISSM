#!/bin/bash
set -eu


# Constants
#
VER="1.6.5"

# Environment
#
export HDF5_DIR="${ISSM_DIR}/externalpackages/petsc/install"
export MPIINC_DIR="${ISSM_DIR}/externalpackages/petsc/install/include"
export NETCDF_DIR="${ISSM_DIR}/externalpackages/netcdf/install"
export PYTHONUSERBASE="${ISSM_DIR}/externalpackages/netcdf-python/install" # This variable and '--user' option supplied to 'setup.py install' are required to install to custom location (source: https://setuptools.readthedocs.io/en/latest/easy_install.html#custom-installation-locations)

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://issm.ess.uci.edu/files/externalpackages/netcdf4-python-${VER}rel.tar.gz" "netcdf4-python-${VER}rel.tar.gz"

# Unpack source
tar -zxvf netcdf4-python-${VER}rel.tar.gz

# Cleanup
rm -rf install src
mkdir install src

# Move source to 'src' directory
mv netcdf4-python-${VER}rel/* src/
rm -rf netcdf4-python-${VER}rel

# Compile and install
cd src
python setup.py build
python setup.py install --user

# Unzip eggs
for egg in $(find "${PYTHONUSERBASE}/lib" -name *.egg); do
	parent_dir=$(dirname ${egg})
	filename=$(basename -- ${egg})
	extension="${filename##*.}"
	filename="${filename%.*}"
	src_dir="${parent_dir}/${filename}"
	unzip ${egg} -d ${src_dir}
done
