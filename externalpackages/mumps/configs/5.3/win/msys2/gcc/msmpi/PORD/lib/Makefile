
# To compile directly, uncomment the line below.
# include ../Make.in
#
# Otherwise, adequate variables for CC, CFLAGS, AR and RANLIB must be passed to 
# make.
#

INCLUDES = -I../include

# NOTE:
#- CFLAGS = $(OPTC) (see ../../src/Makefile; OPTC is set in ../../Makefile.inc)
#
COPTIONS = $(INCLUDES) $(CFLAGS) $(OPTFLAGS)

OBJS = graph.o gbipart.o gbisect.o ddcreate.o ddbisect.o nestdiss.o \
       multisector.o gelim.o bucket.o tree.o \
       symbfac.o interface.o sort.o minpriority.o 

# NOTE: numfac.c read.c mapping.c triangular.c matrix.c kernel.c were not 
# directly used by MUMPS and have been removed from the original SPACE 
# package.

.c.o:
	$(CC) $(COPTIONS) -c $*.c $(OUTC) $*.o

libpord$(LIBEXT):$(OBJS)
	$(CC) -shared $(OUTC) $@ $(OBJS)

clean:
	$(RM) *.o

realclean:
	$(RM) *.o libpord*
