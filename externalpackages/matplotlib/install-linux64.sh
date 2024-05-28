#/bin/bash
set -eu
export GIT_SSL_NO_VERIFY=true 
export CC="gcc -fPIC "
export CXX="g++ -fPIC -L$ISSM_DIR/externalpackages/tcl/install/lib"
export F77="gfortran -fPIC"
export FC="gfortran -fPIC"
export FFLAGS=-ff2c

git clone https://github.com/matplotlib/matplotlib
mv matplotlib src
cd src
python setup.py build 
python setup.py install 
