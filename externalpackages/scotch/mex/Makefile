## Copyright 2004,2007-2009 ENSEIRB, INRIA & CNRS
##
## This file is part of the Scotch software package for static mapping,
## graph partitioning and sparse matrix ordering.
##
## This software is governed by the CeCILL-C license under French law
## and abiding by the rules of distribution of free software. You can
## use, modify and/or redistribute the software under the terms of the
## CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
## URL: "http://www.cecill.info".
## 
## As a counterpart to the access to the source code and rights to copy,
## modify and redistribute granted by the license, users are provided
## only with a limited warranty and the software's author, the holder of
## the economic rights, and the successive licensors have only limited
## liability.
## 
## In this respect, the user's attention is drawn to the risks associated
## with loading, using, modifying and/or developing or reproducing the
## software by the user in light of its specific status of free software,
## that may mean that it is complicated to manipulate, and that also
## therefore means that it is reserved for developers and experienced
## professionals having in-depth computer knowledge. Users are therefore
## encouraged to load and test the software's suitability as regards
## their requirements in conditions enabling the security of their
## systems and/or data to be ensured and, more generally, to use and
## operate it in the same conditions as regards security.
## 
## The fact that you are presently reading this means that you have had
## knowledge of the CeCILL-C license and that you accept its terms.
##

#bindir		= ../../bin
#includedir	= ../../include
#libdir		= ../../lib
#bindir		= ${SCOTCH_DIR}/../src/bin
#includedir	= ${SCOTCH_DIR}/../src/include
#libdir		= ${SCOTCH_DIR}/../src/lib
bindir		= ${SCOTCH_DIR}/bin
includedir	= ${SCOTCH_DIR}/include
libdir		= ${SCOTCH_DIR}/lib

##
##  General inference rules.
##

#include ../Makefile.inc
include ${SCOTCH_DIR}/../src/src/Makefile.inc

%$(OBJ)	:	%.c
				$(CC) $(CFLAGS) -I${SCOTCH_DIR}/../src/src/scotch -I$(includedir) -I${SCOTCH_DIR}/../src/src/libscotch -DSCOTCH_VERSION=\"$(VERSION)\" -c $(<) -o $(@)

%$(MEX)	:	%.c
		 		$(CCM) $(MFLAGS) -I${SCOTCH_DIR}/../src/src/scotch -I$(includedir) -I${SCOTCH_DIR}/../src/src/libscotch -DSCOTCH_VERSION=\"$(VERSION)\" $(<) -o $(@) -L$(libdir) -l$(SCOTCHLIB) -l$(SCOTCHLIB)errexit $(LDFLAGS)

##
##  Project rules.
##

.PHONY				:	mexscotch	mexinstall	clean	realclean

mexscotch			:	clean
					$(MAKE) CFLAGS="$(CFLAGS) -DNOFILEIO" CC="$(CCS)" CCD="$(CCS)" SCOTCHLIB=nfioscotch \
                    Gmapx$(OBJ)
					$(MAKE) CFLAGS="$(CFLAGS) -DNOFILEIO -DMATLAB -I${MATLAB_DIR}/extern/include" CC="$(CCS)" LDFLAGS="$(LDFLAGS) -Wl,-rpath-link,${MATLAB_DIR}/bin/glnxa64 -L${MATLAB_DIR}/bin/glnxa64 -lmex -lmat Gmapx$(OBJ)" SCOTCHLIB=nfioscotch	\
					Gmap$(MEX)

mexinstall			:
					-$(CP) Gmap$(MEX) $(bindir)
					-$(RM) $(bindir)/Gpart$(MEX)
					-$(LN) $(bindir)/Gmap$(MEX) $(bindir)/Gpart$(MEX)

clean				:
					-$(RM) *~ *$(OBJ)
					-$(RM) Gmap$(MEX)

realclean			:	clean

##
##  Todo list.
##

Gmapx$(OBJ)			:	Gmapx.c					\
					${SCOTCH_DIR}/../src/src/libscotch/module.h			\
					${SCOTCH_DIR}/../src/src/libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libnfioscotch$(LIB)		\
					$(libdir)/libnfioscotcherrexit$(LIB)	\
					${SCOTCH_DIR}/../src/src/scotch/gmap.h

Gmap$(MEX)		:	Gmap.c					\
					${SCOTCH_DIR}/../src/src/libscotch/module.h			\
					${SCOTCH_DIR}/../src/src/libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libnfioscotch$(LIB)		\
					$(libdir)/libnfioscotcherrexit$(LIB)	\
					${SCOTCH_DIR}/../src/src/scotch/gmap.h

Gpart$(MEX)		:	Gmap$(MEX)
					-$(RM) Gpart$(MEX)
					-$(LN) Gmap$(MEX) Gpart$(MEX)

