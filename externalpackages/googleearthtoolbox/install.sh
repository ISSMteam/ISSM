#!/bin/bash
set -eu

#Some cleanup
rm -rf install  

#Download code: 
svn checkout http://googleearthtoolbox.googlecode.com/svn/trunk/ install
