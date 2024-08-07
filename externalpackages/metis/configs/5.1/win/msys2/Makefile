# Configuration options.
gdb        = not-set
assert     = not-set
assert2    = not-set
debug      = not-set
gprof      = not-set
openmp     = not-set
prefix     = not-set
gklib_path = not-set
windows    = not-set
shared     = not-set
cc         = not-set
cxx        = not-set
ar         = not-set
gcc-ar     = not-set
ranlib     = not-set
gcc-ranlib = not-set
linker     = not-set
nm         = not-set
objcopy    = not-set
objdump    = not-set
rc         = not-set
strip      = not-set


# Basically proxies everything to the builddir cmake.

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR = build/$(systype)-$(cputype)

# Process configuration options.
CONFIG_FLAGS = -DCMAKE_VERBOSE_MAKEFILE=1
ifeq ($(gklib_path), not-set)
    gklib_path = GKlib
endif
CONFIG_FLAGS += -DGKLIB_PATH=$(abspath $(gklib_path))
ifneq ($(gdb), not-set)
    CONFIG_FLAGS += -DGDB=$(gdb)
endif
ifneq ($(assert), not-set)
    CONFIG_FLAGS += -DASSERT=$(assert)
endif
ifneq ($(assert2), not-set)
    CONFIG_FLAGS += -DASSERT2=$(assert2)
endif
ifneq ($(debug), not-set)
    CONFIG_FLAGS += -DDEBUG=$(debug)
endif
ifneq ($(gprof), not-set)
    CONFIG_FLAGS += -DGPROF=$(gprof)
endif
ifneq ($(openmp), not-set)
    CONFIG_FLAGS += -DOPENMP=$(openmp)
endif
ifneq ($(prefix), not-set)
    CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX=$(prefix)
endif
ifneq ($(windows), not-set)
    CONFIG_FLAGS += -DMSVC=0
endif
ifneq ($(shared), not-set)
    CONFIG_FLAGS += -DSHARED=1 -DCMAKE_SHARED_LINKER_FLAGS=-Wl,--no-undefined
endif
ifneq ($(cc), not-set)
    CONFIG_FLAGS += -DCMAKE_C_COMPILER=$(cc)
endif
ifneq ($(cxx), not-set)
    CONFIG_FLAGS += -DCMAKE_CXX_COMPILER=$(cxx)
endif
ifneq ($(ar), not-set)
    CONFIG_FLAGS += -DCMAKE_AR=$(ar)
endif
ifneq ($(gcc-ar), not-set)
    CONFIG_FLAGS += -DCMAKE_C_COMPILER_AR=$(gcc-ar) -DCMAKE_CXX_COMPILER_AR=$(gcc-ar)
endif
ifneq ($(ranlib), not-set)
    CONFIG_FLAGS += -DCMAKE_RANLIB=$(ranlib)
endif
ifneq ($(gcc-ranlib), not-set)
    CONFIG_FLAGS += -DCMAKE_C_COMPILER_RANLIB=$(gcc-ranlib) -DCMAKE_CXX_COMPILER_RANLIB=$(gcc-ranlib)
endif
ifneq ($(linker), not-set)
    CONFIG_FLAGS += -DCMAKE_LINKER=$(linker)
endif
ifneq ($(nm), not-set)
    CONFIG_FLAGS += -DCMAKE_NM=$(nm)
endif
ifneq ($(objcopy), not-set)
    CONFIG_FLAGS += -DCMAKE_OBJCOPY=$(objcopy)
endif
ifneq ($(objdump), not-set)
    CONFIG_FLAGS += -DCMAKE_OBJDUMP=$(objdump)
endif
ifneq ($(rc), not-set)
    CONFIG_FLAGS += -DCMAKE_RC_COMPILER=$(rc)
endif
ifneq ($(strip), not-set)
    CONFIG_FLAGS += -DCMAKE_STRIP=$(strip)
endif

VERNUM=5.1.0
PKGNAME=metis-$(VERNUM)

define run-config
mkdir -p $(BUILDDIR)
cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all clean install:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more BUILD.txt; \
	else \
	  	make -C $(BUILDDIR) $@ $(MAKEFLAGS); \
	fi

uninstall:
	xargs rm < $(BUILDDIR)/install_manifest.txt

config: distclean
	$(run-config)

distclean:
	rm -rf $(BUILDDIR)

remake:
	find . -name CMakeLists.txt -exec touch {} ';'

dist:
	utils/mkdist.sh $(PKGNAME)

.PHONY: config distclean all clean install uninstall remake dist
