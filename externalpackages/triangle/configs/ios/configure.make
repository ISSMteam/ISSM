CC=${host_triplet}-clang
AR=${host-triplet}-ar
RANLIB=${host-triplet}-ranlib
CSWITCHES = $(CFLAGS) -DNO_TIMER
TRILIBDEFS = -DTRILIBRARY
OBJ_EXT=o
LIB_EXT=a
