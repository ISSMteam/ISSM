#!/bin/bash
set -eu

#Erase symlink
rm -rf install

#symlink to flaim directory
ln -s /home/jschierm/flaim/svn/trunk ./install
