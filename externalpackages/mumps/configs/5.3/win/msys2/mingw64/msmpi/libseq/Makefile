#
#  This file is part of MUMPS 5.3.5, released
#  on Thu Oct 22 09:29:08 UTC 2020
#
all: libmpiseq

.PHONY: all libmpiseq clean

include ../Makefile.inc

libmpiseq: libmpiseq$(PLAT)$(LIBEXT)

# NOTE:
# - We need to use Fortran compiler or otherwise link in Fortran libraries 
#   using C compiler (which does not make a lot of sense, but is possible)
#
libmpiseq$(PLAT)$(LIBEXT): mpi.o mpic.o elapse.o
	$(FC) -shared $(OUTF) $@ mpi.o mpic.o elapse.o
.f.o:
	$(FC) $(OPTF)              -c $*.f $(OUTF)$*.o
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -c $*.c $(OUTC)$*.o

clean:
	$(RM) *.o *$(LIBEXT)
