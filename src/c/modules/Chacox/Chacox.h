/*!\file:  Chacoxx.h
 * \brief header file for Chaco partitioner
 */ 

#ifndef _CHACOX_H
#define _CHACOX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>

#ifdef _HAVE_CHACO_ //only works if dakota library has been compiled in.

#include "chaco.h"

#define    OPT_GLOBAL    0
#define    OPT_LOCAL     1
#define    OPT_VWGTS     2
#define    OPT_EWGTS     3
#define    OPT_ARCH      4
#define    OPT_NDIMS     5
#define    OPT_VMAX      6
#define    OPT_RQI       7
#define    OPT_EIGTOL    8
#define    OPT_SEED      9

#endif

#include "../../classes/classes.h"

/* local prototypes: */
int Chacox( int nvtxs, int *start, int *adjacency,int *vwgts, float *ewgts, float *x, float *y, float *z, short *assignment, double  options[10], int *nparts, double *goal);
int input_parse( char *outassignname, char *outfilename, int *architecture, int *ndims_tot, int mesh_dims[3], int *global_method, int *local_method, 
		int *rqi_flag, int *vmax, int *ndims, int *nprocs, double options[10], int *nparts);
double    chaco_seconds(void);

#undef __FUNCT__ 
#define __FUNCT__  "Chacox"

#endif  /* _CHACOX_H */
