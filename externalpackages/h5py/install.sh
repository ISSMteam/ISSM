#!/bin/bash
set -eu

#needed further along
export HDF5_DIR=$ISSM_DIR/externalpackages/hdf5/install

#Some cleanup
rm -rf install h5py-2.0.1
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/h5py-2.0.1.tar.gz' 'h5py-2.0.1.tar.gz'

#Untar 
tar -zxvf  h5py-2.0.1.tar.gz

#Move h5py to install directory
rm -rf install/*
mv h5py-2.0.1/* install/
rm -rf h5py-2.0.1

#Configure and compile
cd install
python setup.py build –hdf5=$ISSM_DIR/externalpackages/hdf5/install
python setup.py install
