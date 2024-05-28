#!/bin/bash
set -eu


## Constants
#
PREFIX="${ISSM_DIR}/externalpackages/semic/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX}/lib

# Download source
git clone https://github.com/mkrapp/semic.git src

if which ifort >/dev/null; then
	FC="ifort"
	FFLAGS="-traceback -check all" # -O2 is default 
else
	FC="gfortran"
	if [ `uname` == "Darwin" ]; then
		FFLAGS="-fcheck=all -arch $(uname -m)"
	else
		FFLAGS=""
	fi
fi

# Compile and install semic module utils.f90
cd src/
(
cat << EOF
LIB_EXT=a
FC=$FC
FFLAGS=$FFLAGS
install: libutils.\$(LIB_EXT)
	cp libutils.\$(LIB_EXT) ${PREFIX}/lib
	cp utils.mod ${PREFIX}
OBJECTS= utils.o
libutils.\$(LIB_EXT): \$(OBJECTS)
	ar -r libutils.\$(LIB_EXT) \$(OBJECTS) 
	ranlib libutils.\$(LIB_EXT) 
%.o: %.f90
	\$(FC) \$(FFLAGS) -fPIC -c $< -o \$@
clean: 
	rm -rf *.o *.\$(LIB_EXT)
EOF
) > Makefile
make

# Apply patch surface_physics
patch surface_physics.f90 < ../configs/surface_physics.f90.patch

# Compile semic module surface_physics.f90
(
cat << EOF
LIB_EXT=a
FC=$FC
FFLAGS=$FFLAGS
install: libsurface_physics.\$(LIB_EXT)
	cp libsurface_physics.\$(LIB_EXT) ${PREFIX}/lib
	cp surface_physics.mod ${PREFIX}
OBJECTS= surface_physics.o
libsurface_physics.\$(LIB_EXT): \$(OBJECTS)
	ar -r libsurface_physics.\$(LIB_EXT) \$(OBJECTS) 
	ranlib libsurface_physics.\$(LIB_EXT) 
%.o: %.f90
	\$(FC) \$(FFLAGS) -fPIC -c $< -o \$@
clean: 
	rm -rf *.o *.\$(LIB_EXT)
EOF
) > Makefile
make
