#!/bin/bash
set -eu

#Some cleanup
rm -rf install src cccl-0.03
mkdir install 
mkdir install/bin

#Move cccl into install directory
cp issm/cccl install/bin
chmod 755 install/bin/cccl
