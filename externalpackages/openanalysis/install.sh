#!/bin/bash
set -eu

#Some cleanup
rm -rf openanalysis

#download openanalysis
svn co http://svn.berlios.de/svnroot/repos/openanalysis/OpenAnalysis/trunk openanalysis

#Configure
cd openanalysis

make -f Makefile.quick all
make -f Makefile.quick install
