#!/bin/bash
set -eu

#Erase install
rm -rf install  src 
mkdir  install

#Download
git clone git@github.com:danmar/cppcheck.git src

#compile
cd src
if [ $# -eq 0 ]; then
	make MATCHCOMPILER=yes FILESDIR=$ISSM_DIR/externalpackages/cppcheck/install HAVE_RULES=yes CXXFLAGS="-O2 -DNDEBUG -Wall -Wno-sign-compare -Wno-unused-function"
else 
	make -j $1 MATCHCOMPILER=yes FILESDIR="$ISSM_DIR/externalpackages/cppcheck/install" HAVE_RULES=yes CXXFLAGS="-O2 -DNDEBUG -Wall -Wno-sign-compare -Wno-unused-function"
fi  
cd ..
