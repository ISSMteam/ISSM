#!/bin/bash
set -eu

export CC="gcc -fPIC"
export CXX="g++ -fPIC"
export F77="gfortran -fPIC"
export FC="gfortran -fPIC"
export FFLAGS=-ff2c

#clean up
rm -rf scipy
#download scipy
export GIT_SSL_NO_VERIFY=true 
git clone https://github.com/scipy/scipy.git

#install scipy
cd scipy
export  BLAS_SRC=$ISSM_DIR/externalpackages/blas/install/lib
export  BLAS=$ISSM_DIR/externalpackages/blas/install/lib
export  LAPACK_SRC=$ISSM_DIR/externalpackages/lapack/install/lib
export  LAPACK=$ISSM_DIR/externalpackages/lapack/install/lib

python setup.py build
python setup.py install
cd ..
python -c "import scipy; print 'Installed SciPy', scipy.__version__"
#python -c "import scipy; scipy.test()"
