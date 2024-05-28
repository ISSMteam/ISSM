#!/bin/bash
set -eu

# On OSX 10.6, fgortran gets installed in /usr/local/gfortran 
#export PATH="/usr/local/gfortran/bin/:$PATH"
export F77="/usr/local/gfortran/bin/x86_64-apple-darwin10-gfortran"
export CC="/usr/bin/gcc"
export CXX="/usr/bin/g++"

#download scipy
git clone https://github.com/scipy/scipy.git

#install scipy
cd scipy
python setup.py build
python setup.py install
cd ..
python -c "import scipy; print 'Installed SciPy', scipy.__version__"
#python -c "import scipy; scipy.test()"
