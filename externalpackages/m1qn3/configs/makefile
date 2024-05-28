include ./configure.make

ifeq "$(origin AR)" "undefined"
	AR = ar
endif

ifeq "$(origin RANLIB)" "undefined"
	RANLIB = ranlib
endif

OBJECTS=${LIBNAME}.o
SOURCES=${LIBNAME}.f

objects: $(OBJECTS)

shared: lib${LIBNAME}.$(SHARED_LIB_EXT)

static: lib${LIBNAME}.$(STATIC_LIB_EXT)

${LIBNAME}.o: $(SOURCES)
	${FC} ${FFLAGS} -fpic -c $< -o $@

lib${LIBNAME}.a: $(OBJECTS)
	$(AR) -cr $@ $(OBJECTS)
	$(RANLIB) $@

lib${LIBNAME}.dll: $(OBJECTS)
	${FC} ${FFLAGS} -shared -c $(SOURCES) -o $@

lib${LIBNAME}.dylib: $(OBJECTS)
	${FC} ${FFLAGS} -dynamiclib -install_name ${PREFIX}/lib/$@ -c $(SOURCES) -o $@

lib${LIBNAME}.so: $(OBJECTS)
	${FC} ${FFLAGS} -shared -c $(SOURCES) -o $@

clean: 
	rm -rf *.o *.$(SHARED_LIB_EXT) *.$(STATIC_LIB_EXT)
