#!/bin/bash
set -eu

#Note of caution:  stop after boostrap phase, and run 
#bjam --debug-configuration, to figure out which paths boost is using to include 
#python. make sure everyone of these paths is covered by python. If not, just make 
#symlinks in externalpackages/python to what boost is expecting. Ther is NO WAY 
#to get the boost library to include python support without doing that. 

#Some cleanup
rm -rf install boost_1_55_0 src
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/boost_1_55_0.tar.gz' 'boost_1_55_0.tar.gz'

#Untar 
tar -zxvf  boost_1_55_0.tar.gz

#Move boost into install directory
mv boost_1_55_0/* src
rm -rf boost_1_55_0

patch src/boost/atomic/detail/cas128strong.hpp ./configs/1.55/cas128strong.hpp.patch
patch src/boost/atomic/detail/gcc-atomic.hpp ./configs/1.55/gcc-atomic.hpp.patch
patch src/tools/build/v2/user-config.jam ./configs/1.55/user-config.jam.patch
patch src/tools/build/v2/tools/darwin.jam ./configs/1.55/darwin.jam.patch
patch src/tools/build/v2/tools/darwin.py ./configs/1.55/darwin.py.patch

#Configure and compile
cd src 
./bootstrap.sh \
	--prefix="$ISSM_DIR/externalpackages/boost/install" \
	--with-python=python 

#Compile boost
# Need gcc with iconv installed in a location that has been added to your path
#./bjam toolset=darwin-std0x link=static install

export CC=/usr/local/gfortan/bin/gcc 
export CXX=/usr/local/gfortran/bin/g++
./b2 toolset=clang cxxflags=-stdlib=libstdc++ linkflags=-stdlib=libstdc++ -j2 architecture=x86 variant=release link=static threading=multi install

#put bjam into install also: 
mkdir ../install/bin
cp bjam ../install/bin
