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

CSWITCHES=$(CFLAGS) -DNO_TIMER -fpic
OBJECTS=triangle.o
SOURCES=triangle.c triangle.h
TRILIBDEFS=-DTRILIBRARY

all: libtriangle.$(STATIC_LIB_EXT) libtriangle.$(SHARED_LIB_EXT)

objects: $(OBJECTS)

shared: libtriangle.$(SHARED_LIB_EXT)

static: libtriangle.$(STATIC_LIB_EXT)

triangle.o: $(SOURCES)
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -c triangle.c

libtriangle.a: $(OBJECTS)
	$(AR) -cr $@ $(OBJECTS)
	$(RANLIB) $@

libtriangle.dll: $(SOURCES)
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -shared -o $@ triangle.c

libtriangle.dylib: $(OBJECTS)
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -dynamiclib -install_name ${PREFIX}/lib/$@ -o $@ triangle.c

libtriangle.lib: $(OBJECTS)
	lib -out:libtriangle.$(STATIC_LIB_EXT) $(OBJECTS)

libtriangle.so: $(SOURCES)
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -shared -o $@ triangle.c

clean:
	rm -rf *.$(LIB_EXT) *.o *.LIB *.LST *.obj *.BAK
