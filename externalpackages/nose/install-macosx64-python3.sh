#!/bin/bash
#Install Python nose module

rm -rf src  install

svn checkout http://python-nose.googlecode.com/svn/branches/py3k
mv py3k src

cd src
python ./setup.py build
python ./setup.py install

#to be flagged by jenkins, we create an empty install dir: 
cd ../
mkdir install
touch install/emptyfile
