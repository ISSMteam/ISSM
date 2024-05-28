#!/bin/bash
set -eu

#clean up
rm -rf cython

#download cython first
git clone https://github.com/cython/cython.git

#install cython
cd cython
python setup.py install
