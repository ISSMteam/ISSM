#Adaptation of Triangle makefile to compile only the triangle.c file, 
#into a standalone library that can then be used to compile mex modules

include ./configure.make

ifeq "$(origin CC)" "undefined"
	CC = cc
endif

ifeq "$(origin AR)" "undefined"
	AR = ar
endif

ifeq "$(origin RANLIB)" "undefined"
	RANLIB = ranlib
endif

TRILIBDEFS = -DTRILIBRARY

all: triangle.$(LIB_EXT)

OBJECTS=triangle.o

triangle.$(LIB_EXT): $(OBJECTS)
	$(AR) cr  triangle.$(LIB_EXT) $(OBJECTS)
	$(RANLIB) triangle.$(LIB_EXT)

triangle.o: triangle.c triangle.h
	$(CXX) $(CXXSWITCHES) $(TRILIBDEFS) -c triangle.c

clean: 
	rm -rf *.$(LIB_EXT) *.o *.LIB *.LST *.obj *.BAK
