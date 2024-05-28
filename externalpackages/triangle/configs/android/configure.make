# This makefile configures build process to cross-compile to the android platform.
# The binary tools referenced below are specifically configured to target armeabi-v7a.
# Furthermore, the compilers (which are simply wrappers around GNU GCC) are set to
# produce binaries that are EABI compliant.
#
# Note that the AAPCS standard defines 'EABI' as a moniker used to specify
# a _family_ of similar but distinct ABIs. Android follows the little-endian
# ARM GNU/Linux ABI as documented in the following document:
#
# http://www.codesourcery.com/gnu_toolchains/arm/arm_gnu_linux_abi.pdf

CC=${host_triplet}-gcc
AR=${host_triplet}-ar
RANLIB=${host_triplet}-ranlib
CSWITCHES = $(CFLAGS) -DNO_TIMER 
TRILIBDEFS = -DTRILIBRARY
OBJ_EXT=o
LIB_EXT=a
