#!/bin/bash
set -eu


# TODO:
# - Introduce build step to $ISSM_DIR/jenkins/jenkins.sh to compile Fortran code in $ISSM_DIR/src/c/modules/GiaDefelectionCorex/ to C with f2c
#       - Then, revert $ISSM_DIR/externalpackages/emscripten/install.sh to r24306 and test clean build
#       - When builtin support for Fortran is available, remove build step
#

## Constants
#
VER="latest" # Set this to "latest", or last tag that works in case of failure

## Environment
#
PREFIX="${ISSM_DIR}/externalpackages/emscripten/install"

# Cleanup
rm -rf ${PREFIX}

# Get the emsdk repo
git clone https://github.com/emscripten-core/emsdk.git

# Create $PREFIX directory
mkdir -p ${PREFIX}

# Move source to $PREFIX directory
mv emsdk/* ${PREFIX}
rm -rf emsdk

cd ${PREFIX}

# Download and install the latest SDK tools.
./emsdk install ${VER}

# Make the "latest" SDK "active" for the current user. (writes ~/.emscripten
# file)
./emsdk activate ${VER}

# Activate PATH and other environment variables in the current terminal
source ./emsdk_env.sh
