/*!\file: chaco.h
 * \brief prototypes for chaco library
 */ 

#ifndef _CHACO_INCLUDES_H_
#define  _CHACO_INCLUDES_H_

#include "./defs.h"
#include "./params.h"

extern "C" int       interface( int       nvtxs,		/* number of vertices in full graph */
		int      *start,		/* start of edge list for each vertex */
		int      *adjacency,		/* edge list data */
		int      *vwgts,		/* weights for all vertices */
		float    *ewgts,		/* weights for all edges */
		float    *x, 
		float    *y, 
		float    *z,		/* coordinates for inertial method */
		char     *outassignname,	/* name of assignment output file */
		char     *outfilename,		/* output file name */
		short    *assignment,		/* set number of each vtx (length n) */
		int       architecture,		/* 0 => hypercube, d => d-dimensional mesh */
		int       ndims_tot,		/* total number of cube dimensions to divide */
		int       mesh_dims[3],		/* dimensions of mesh of processors */
		double   *goal,			/* desired set sizes for each set */
		int       global_method,	/* global partitioning algorithm */
		int       local_method,		/* local partitioning algorithm */
		int       rqi_flag,		/* should I use RQI/Symmlq eigensolver? */
		int       vmax,			/* how many vertices to coarsen down to? */
		int       ndims,		/* number of eigenvectors (2^d sets) */
		double    eigtol,		/* tolerance on eigenvectors */
		long      seed);			/* for random graph mutations */

extern "C" void      read_params(FILE* pfile);	/* file with new user parameters */
extern "C" void      smalloc_stats(void);
extern "C" double*   smalloc(unsigned int n); /* number of bytes to be allocated */

#endif //ifndef _CHACO_INCLUDES_H_

