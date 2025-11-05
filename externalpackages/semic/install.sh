#!/bin/bash
set -eu


## Constants
#
VER="1.2"

PREFIX="${ISSM_DIR}/externalpackages/semic/install" # Set to location where external package should be installed

# Cleanup
rm -rf ${PREFIX} src
mkdir -p ${PREFIX} ${PREFIX}/lib src

# Download source
$ISSM_DIR/scripts/DownloadExternalPackage.sh "https://github.com/mkrapp/semic/archive/refs/tags/v${VER}.tar.gz" "semic-${VER}.tar.gz"

# Unpack source
tar -xvzf semic-${VER}.tar.gz

# Move source to 'src' directory
mv semic-${VER}/* src
rm -rf semic-${VER}

# Compile and install semic module utils.f90
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
