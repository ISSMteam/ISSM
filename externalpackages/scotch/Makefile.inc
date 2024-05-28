EXE	=
LIB	= .a
OBJ	= .o
MEX	= .mexa64

MAKE	= make
AR	= ar
ARFLAGS	= -ruv
CAT	= cat
CCS	= gcc
CCP     = mpicc
CCD     = mpicc
CCM	= ${MATLAB_DIR}/bin/mex
#CFLAGS	= -m64 -O3 -std=c99 -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Du_int32_t=uint32_t -Du_int64_t=uint64_t -DCOMMON_TIMING_OLD -DMATLAB -fPIC -I${MATLAB_DIR}/extern/include
CFLAGS	= -m64 -O3 -std=c99 -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Du_int32_t=uint32_t -Du_int64_t=uint64_t -DCOMMON_TIMING_OLD -fPIC
#LDFLAGS	= -lz -lm -lrt -L${MATLAB_DIR}/bin/glnxa64 -lmex -lmat
#LDFLAGS	= -lz -lm -lrt -Wl,-rpath-link,${MATLAB_DIR}/bin/glnxa64 -L${MATLAB_DIR}/bin/glnxa64 -lmex -lmat
LDFLAGS	= -lz -lm -lrt
MFLAGS	= -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Du_int32_t=uint32_t -Du_int64_t=uint64_t -DCOMMON_TIMING_OLD -DNOFILEIO -DMATLAB -I${MATLAB_DIR}/extern/include -largeArrayDims
CP	= cp
LEX	= flex
LN	= ln
MKDIR	= mkdir
MV	= mv
RANLIB	= ranlib
YACC	= yacc
