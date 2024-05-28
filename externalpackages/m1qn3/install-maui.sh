#!/bin/bash
set -eu

#Some cleanup 
rm -rf install src m1qn3-3.3-distrib
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/m1qn3-3.3-distrib.tgz' 'm1qn3-3.3-distrib.tgz'

#Untar 
tar -xzf m1qn3-3.3-distrib.tgz
mv m1qn3-3.3-distrib src

FC="ftn"

#Compile m1qn3
cd src/src/
(
cat << EOF
LIB_EXT=a
FC=$FC
install: libm1qn3.\$(LIB_EXT)
	cp libm1qn3.\$(LIB_EXT) ../../install/
OBJECTS= m1qn3.o
libm1qn3.\$(LIB_EXT): \$(OBJECTS)
	ar -r libm1qn3.\$(LIB_EXT) \$(OBJECTS) 
	ranlib libm1qn3.\$(LIB_EXT) 
%.o: %.f
	\$(FC) \$(FFLAGS) -fPIC -c $< -o \$@
clean: 
	rm -rf *.o *.\$(LIB_EXT)
EOF
) > Makefile
make

#compile ddot
cd ../blas
(
cat << EOF
LIB_EXT=a
FC=$FC
install: libddot.\$(LIB_EXT)
	cp libddot.\$(LIB_EXT) ../../install/
OBJECTS= ddot.o
libddot.\$(LIB_EXT): \$(OBJECTS)
	ar -r libddot.\$(LIB_EXT) \$(OBJECTS) 
	ranlib libddot.\$(LIB_EXT) 
%.o: %.f
	\$(FC) \$(FFLAGS) -fPIC -c $< -o \$@
clean: 
	rm -rf *.o *.\$(LIB_EXT)
EOF
) > Makefile
make
