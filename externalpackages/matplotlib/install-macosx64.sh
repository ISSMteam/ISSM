#/bin/bash
set -eu

rm -rf src deps 

export C_INCLUDE_PATH="$C_INCLUDE_PATH:/usr/include/c++/4.2.1"
export GIT_SSL_NO_VERIFY=true 
export CC="gcc -fPIC -arch x86_64 -I/usr/include/c++/4.2.1"
export CXX="g++ -fPIC -arch x86_64 -I/usr/include/c++/4.2.1"
export F77="gfortran -fPIC -arch x86_64 -I/usr/include/c++/4.2.1"
export FC="gfortran -fPIC -arch x86_64 -I/usr/include/c++/4.2.1"
export FFLAGS="-ff2c -arch x86_64 -I/usr/include/c++/4.2.1"
export ARCHFLAGS="-arch x86_64 -I/usr/include/c++/4.2.1"

git clone https://github.com/matplotlib/matplotlib
mv matplotlib src
mkdir deps
cd src

#only try this if the classic python setup.py build approach does not work. The approach in the next 3 lines can 
#trigger issues of permissions with freetype. Sometimes it is better to independently install freetype2 (from the 
#issm externalpackages) as sudo (or root), so that the python script can detect its existence automatically.
#sudo make -f make.osx PREFIX=$ISSM_DIR/externalpackages/matplotlib/deps PYVERSION=$pythonversion fetch deps mpl_install_std
#sudo make -f make.osx PREFIX=$ISSM_DIR/externalpackages/matplotlib/deps PYVERSION=$pythonversion mpl_install_std
#python -c "import matplotlib; print 'Installed matplotlib', matplotlib.__version__, matplotlib.__file__"

#comments: try exporting this first before calling python setup.py build
CFLAGS=-mmacosx-version-min=10.7 

#to be tried:  first get freetype and zlib and libpng installed in sudo mode
python setup.py build
python setup.py install
