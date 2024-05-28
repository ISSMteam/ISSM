#!/bin/bash
set -eu

#Some cleanup
rm -rf angel

#download 
svn co -r 82 https://angellib.svn.sourceforge.net/svnroot/angellib/trunk angel

#Compile
cd angel 
make
