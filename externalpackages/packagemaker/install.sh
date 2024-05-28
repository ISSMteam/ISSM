#!/bin/bash
set -eu

#Erase symlink
rm -rf install

#Select or create a new simlink
ln -s /Developer/Applications/Utilities/PackageMaker.app/Contents/MacOS/ install
