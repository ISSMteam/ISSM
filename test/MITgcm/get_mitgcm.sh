#!/bin/bash

if [ -e ~/.bashrc ]; then 
	source ~/.bashrc
fi

# Download fresh copy of MITgcm
cd ../MITgcm/
if [ ! -d install ]; then
    source install.sh
else
    cd install
    git pull
    cd ..
fi
