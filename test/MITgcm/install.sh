#!/bin/bash


# Cleanup
rm -rf build/* checkpoint* install

################################################################################

# # Download source
# git clone --depth=1 https://github.com/MITgcm/MITgcm.git

# # Move source to 'install' directory
# mv MITgcm install

################################################################################

# Comment out the above, uncomment the following and set to a specific tagged release if the MITgcm repo head is buggy
#

# Constants
#
#VER="68h"
VER="68i"

# Download source
wget https://github.com/MITgcm/MITgcm/archive/refs/tags/checkpoint${VER}.tar.gz

# Uncompress source
tar -xvzf checkpoint${VER}.tar.gz

# Move source to 'install' directory
mv MITgcm-checkpoint${VER} install
