#!/bin/bash
set -eu

#Some cleanup
rm -rf xaifBooster

#download 
svn co -r 125  http://hpc.svn.rice.edu/r/xaifBooster/trunk xaifBooster

#Compile xaifBooster
cd xaifBooster
make
