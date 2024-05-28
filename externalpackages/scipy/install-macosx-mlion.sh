#!/bin/bash
set -eu

export CC=gcc
export CXX=g++
export FFLAGS=-ff2c
	
#download scipy
git clone https://github.com/scipy/scipy.git

#install scipy
cd scipy
export  BLAS_SRC=$ISSM_DIR/externalpackages/blas/install/lib
export  BLAS=$ISSM_DIR/externalpackages/blas/install/lib
export  LAPACK_SRC=$ISSM_DIR/externalpackages/lapack/install/lib
export  LAPACK=$ISSM_DIR/externalpackages/lapack/install/lib

#install scipy
python setup.py build
python setup.py install
cd ..
python -c "import scipy; print 'Installed SciPy', scipy.__version__"
#python -c "import scipy; scipy.test()"
