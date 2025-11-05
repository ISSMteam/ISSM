#!/bin/bash

./configure \
    --prefix=$ISSM_DIR \
    --with-matlab-dir="/Applications/MATLAB_R2019b.app/" \
    --with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
    --with-mpi-include="$ISSM_DIR/externalpackages/mpich/install/include" \
    --with-mpi-libflags="-L$ISSM_DIR/externalpackages/mpich/install/lib/ -lmpich" \
    --with-petsc-dir="$ISSM_DIR/externalpackages/petsc/install" \
    --with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
    --with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install/" \
    --with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install/" \
    --with-m1qn3-dir="$ISSM_DIR/externalpackages/m1qn3/install" \
    --with-ocean="yes" \
    --with-numthreads=8 \
    --with-cxxoptflags=" -D_DO_NOT_LOAD_GLOBALS_ -g -O2 -std=c++11" \
    --enable-development \
    --enable-debugging
