#!/bin/bash
set -eu

#AppScan install directory. Just symlink to your existing AppScan software

#Erase symlink
rm -rf install

#Select or create a new simlink
ln -s /opt/IBM/AppScan_Source ./install
