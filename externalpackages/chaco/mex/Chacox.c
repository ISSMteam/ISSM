/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "params.h"

#define THISFUNCTION "Chacox"

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

int      input_parse(
	char     *outassignname,	/* name of assignment output file */
	char     *outfilename,		/* name of file for outputing run results */
	int      *architecture,		/* 0=> hypercube, d=> d-dimensional mesh */
	int      *ndims_tot,		/* target number of hypercube dimensions */
	int       mesh_dims[3],		/* mesh dimensions */
	int      *global_method,	/* what global partitioning strategy to use? */
	int      *local_method,		/* what local refinement strategy to use? */
	int      *rqi_flag,		/* should I use multilevel eigensolver? */
	int      *vmax,			/* if so, how far should I coarsen? */
	int      *ndims,		/* number of divisions at each stage */
	int      *nprocs,		/* number of processors being divided into */
	double   options[10],	/* architecture and partitioning options */
	int      *nparts		/* number of parts options */
);


int Chacox(
	int       nvtxs,		/* number of vertices in graph */
	int      *start,		/* start of edge list for each vertex */
	int      *adjacency,	/* edge list data */
	int      *vwgts,		/* weights for all vertices */
	float    *ewgts,		/* weights for all edges */
	float    *x,
	float    *y,
	float    *z,			/* coordinates for inertial method */
	short    *assignment,	/* set number of each vtx (length nvtxs+1) */
	double   options[10],	/* architecture and partitioning options */
	int      *nparts,		/* number of parts options */
	double   *goal			/* desired set sizes */
)
{
	extern int Using_Main;	/* is main routine being called? */
	extern char *PARAMS_FILENAME;	/* name of file with parameter updates */
	extern double EIGEN_TOLERANCE;	/* tolerance for eigen calculations */
	extern int OUTPUT_ASSIGN;	/* whether to write assignment to file */
	extern int DEBUG_MEMORY;	/* debug memory allocation and freeing? */
	extern int DEBUG_TRACE;	/* trace main execution path */
	extern int DEBUG_PARAMS;	/* debug flag for reading parameters */
	extern long RANDOM_SEED;	/* seed for random number generators */
	extern int ECHO;		/* controls amount of output */
	extern int PROMPT;		/* prompt for input or not? */
	extern int PRINT_HEADERS;	/* print lines for output sections? */
	extern int MATCH_TYPE;      /* matching routine to call */
	extern double input_time;	/* times data file input */
	extern double start_time;	/* time partitioning starts */
	FILE     *params_file;	/* file with parameter value updates */
	int       global_method;	/* global partitioning method */
	int       local_method;	/* local partitioning method */
	double    eigtol;		/* tolerance in eigenvector calculation */
	int       ndims;		/* dimension of recursive partitioning */
	int       architecture;	/* 0 => hypercube, d => d-dimensional mesh */
	int       ndims_tot;	/* total number of cube dimensions to divide */
	int       mesh_dims[3];	/* dimensions of mesh of processors */
	long      seed;		/* for random graph mutations */
	int       rqi_flag;		/* use RQI/Symmlq eigensolver? */
	int       vmax;		/* if so, how many vertices to coarsen down to? */
	char      outassignname[NAME_LENGTH];	/* assignment output file name */
	char      outfilename[NAME_LENGTH];	/* name of output file */
	char     *outassignptr;	/* name or null pointer for output assignment */
	char     *outfileptr;	/* name or null pointer for output file */
	int       nprocs;		/* number of processors being divided into */
	double    time;		/* timing marker */
	int       flag;		/* return code from input routines */
	double   *smalloc();	/* safe version of malloc */
	double    seconds();	/* returns elapsed time in seconds */
	int       sfree(), interface(), affirm();
	void      input_queries(), smalloc_stats(), read_params(), clear_timing();

	int i,tvwgt;
	double tgoal;

	if (DEBUG_TRACE > 0) {
		printf("<Entering main>\n");
	}

	if (PRINT_HEADERS) {
		printf("\n                    Chaco 2.0\n");
		printf("          Sandia National Laboratories\n\n");
	}

	Using_Main = TRUE;
	params_file = fopen(PARAMS_FILENAME, "r");
	if (params_file == NULL && DEBUG_PARAMS > 1) {
		printf("Parameter file `%s' not found; using default parameters.\n",
			   PARAMS_FILENAME);
	}

	start_time = time = seconds();

	read_params(params_file);

	flag = input_parse(outassignname, outfilename,
			  &architecture, &ndims_tot, mesh_dims,
			  &global_method, &local_method, &rqi_flag, &vmax, &ndims, &nprocs,
			  options, nparts);
	if (flag)
		return(flag);

	if (OUTPUT_ASSIGN > 0)
		outassignptr = outassignname;
	else
		outassignptr = NULL;

	if (ECHO < 0)
		outfileptr = outfilename;
	else
		outfileptr = NULL;

	if ((int)options[OPT_VWGTS] && vwgts) {
		printf("%s -- Applying weights for %d vertices.\n",
			   THISFUNCTION,nvtxs);
		tvwgt = 0.;
		for (i=0; i<nvtxs; i++)
			tvwgt += vwgts[i];
	}
	else {
		tvwgt = nvtxs;
		if      ( (int)options[OPT_VWGTS] && !vwgts)
			printf("%s -- Vertex weight flag=%d, but no vertex weights specified.\n",
				   THISFUNCTION,options[OPT_VWGTS]);
		else if (!(int)options[OPT_VWGTS] &&  vwgts)
			printf("%s -- Vertex weight flag=%d, so specified vertex weights ignored.\n",
				   THISFUNCTION,options[OPT_VWGTS]);
	}

	if ((int)options[OPT_EWGTS] && ewgts) {
		printf("%s -- Applying weights for %d edges.\n",
			   THISFUNCTION,start[nvtxs]/2);
	}
	else {
		if      ( (int)options[OPT_EWGTS] && !ewgts)
			printf("%s -- Edge weight flag=%d, but no edge weights specified.\n",
				   THISFUNCTION,options[OPT_EWGTS]);
		else if (!(int)options[OPT_EWGTS] &&  ewgts)
			printf("%s -- Edge weight flag=%d, so specified edge weights ignored.\n",
				   THISFUNCTION,options[OPT_EWGTS]);
	}

    if (goal) {
        printf("%s -- Applying goals for %d sets.\n",
               THISFUNCTION,nprocs);
        tgoal = 0.;
        for (i=0; i<nprocs; i++)
            tgoal += goal[i];
        for (i=0; i<nprocs; i++)
            goal[i] *= (double)tvwgt/tgoal;
    }

	input_time += seconds() - time;

	if (options[OPT_EIGTOL] > 0)
		eigtol = options[OPT_EIGTOL];
	else
		eigtol = EIGEN_TOLERANCE;
	if ((int)options[OPT_SEED] > 0)
		seed = (int)options[OPT_SEED];
	else
		seed = RANDOM_SEED;

/*  Chaco numbers vertices from 1 and the Matlab sparse data structure
	numbers rows from 0, so increment the row indices for each column. */

	for (i=0; i<start[nvtxs]; adjacency[i++]++);

	printf("\n%s -- Calling Chaco interface:\n\n",
		   THISFUNCTION);
	flag = interface(nvtxs, start, adjacency,
		  ((int)options[OPT_VWGTS] && vwgts ? vwgts : NULL),
		  ((int)options[OPT_EWGTS] && ewgts ? ewgts : NULL),
		  x, y, z,
		  outassignptr, outfileptr,
		  assignment,
		  architecture, ndims_tot, mesh_dims, goal,
		  global_method, local_method, rqi_flag, vmax, ndims,
		  eigtol, seed);
	printf("\n%s -- Chaco interface returning flag=%d.\n",
		   THISFUNCTION,flag);

/*  Reset adjacency matrix in case calling function needs it.  */

	for (i=0; i<start[nvtxs]; adjacency[i++]--);

	if (DEBUG_MEMORY > 0) {
		printf("\n");
		smalloc_stats();
	}

	if (params_file != NULL)
		fclose(params_file);

	if (DEBUG_TRACE > 1) {
		printf("<Leaving main>\n");
	}
	
	return(0);
}


int      input_parse(
	char     *outassignname,	/* name of assignment output file */
	char     *outfilename,		/* name of file for outputing run results */
	int      *architecture,		/* 0=> hypercube, d=> d-dimensional mesh */
	int      *ndims_tot,		/* target number of hypercube dimensions */
	int       mesh_dims[3],		/* mesh dimensions */
	int      *global_method,	/* what global partitioning strategy to use? */
	int      *local_method,		/* what local refinement strategy to use? */
	int      *rqi_flag,		/* should I use multilevel eigensolver? */
	int      *vmax,			/* if so, how far should I coarsen? */
	int      *ndims,		/* number of divisions at each stage */
	int      *nprocs,		/* number of processors being divided into */
	double   options[10],	/* architecture and partitioning options */
	int      *nparts		/* number of parts options */
)
{
	extern int SEQUENCE;	/* sequence instead of partition graph? */
	extern int ARCHITECTURE;	/* 0=> hypercube, d=> d-dimensional mesh */
	extern int OUTPUT_ASSIGN;	/* write assignments to file? */
	extern int ECHO;		/* copy input to screen? results to file? */
	extern int DEBUG_TRACE;	/* trace main execution path */
	extern int PROMPT;		/* prompt for input? */
	extern int MATCH_TYPE;      /* max-matching routine to call */
	int       eigensolver;	/* which kind of eigensolver to use */

	if (DEBUG_TRACE > 0) {
		printf("<Entering input_parse>\n");
	}

	if (PROMPT) {
		printf("Parallel machine architecture:\n");
		printf("  (0) Hypercube\n");
		printf("  (1) One-dimensional mesh\n");
		printf("  (2) Two-dimensional mesh\n");
		printf("  (3) Three-dimensional mesh\n");
	}
	*architecture = (int)options[OPT_ARCH];
	if (*architecture < 0 || *architecture > 3) {
		printf("%s -- Architecture %d must be between 0 and 3.\n",
			   THISFUNCTION,options[OPT_ARCH]);
		return(-1);
	}

	/* Name output assignment file. */
	if (PROMPT)
		printf("Assignment output file: ");
	outassignname = NULL;

	/* Name output results file. */
	if (PROMPT)
		printf("File name for saving run results: ");
	outfilename = NULL;

	/* Initialize the method flags */
	*rqi_flag = 0;
	*global_method = 0;

	/* Get global method, if any. */
	if (SEQUENCE) {
		*global_method = 2;
	}
	else {
		if (PROMPT) {
			printf("Global partitioning method:\n");
			printf("  (1) Multilevel-KL\n");
			printf("  (2) Spectral\n");
			printf("  (3) Inertial\n");
			printf("  (4) Linear\n");
			printf("  (5) Random\n");
			printf("  (6) Scattered\n");
			printf("  (7) Read-from-file\n");
		}
		*global_method = (int)options[OPT_GLOBAL];
		if (*global_method < 1 || *global_method > 7) {
			printf("%s -- Global method %d must be between 1 and 7.\n",
				   THISFUNCTION,options[OPT_GLOBAL]);
			return(-1);
		}
	}

	if (*global_method == 7) {	/* Name and open input assignment file. */
		if (PROMPT)
			printf("Assignment input file: ");
	}

	else if (*global_method == 3) {
		if (PROMPT)
			printf("Geometry input file name: ");
	}

	else if (*global_method == 2) {
		if (PROMPT) {
			printf("Eigensolver:\n");
			printf("  (1) Multilevel RQI/Symmlq\n");
			printf("  (2) Lanczos\n"); 
		}
		eigensolver = (int)options[OPT_RQI];
		if (eigensolver < 0 || eigensolver > 2) {
			printf("%s -- RQI/Symmlq flag %d must be between 0 and 2.\n",
				   THISFUNCTION,options[OPT_RQI]);
			return(-1);
		}
		if (eigensolver == 1) {
			if (MATCH_TYPE == 5) {	/* geometric matching */
				if (PROMPT)
					printf("Geometry input file name: ");
			}
			*rqi_flag = 1;
			if (PROMPT)
				printf("Number of vertices to coarsen down to: ");
			*vmax = (int)options[OPT_VMAX];
			if (*vmax <= 0) {
				printf("%s -- Vmax %d must be greater then 0.\n",
					   THISFUNCTION,options[OPT_VMAX]);
				return(-1);
			}
		}
		else if (eigensolver == 0 || eigensolver == 2) {
			*rqi_flag = 0;
		}
	}

	else if (*global_method == 1) {
		if (MATCH_TYPE == 5) {		/* geometric matching */
			if (PROMPT)
				printf("Geometry input file name: ");
		}
		if (PROMPT)
			printf("Number of vertices to coarsen down to: ");
		*vmax = (int)options[OPT_VMAX];
		if (*vmax <= 0) {
			printf("%s -- Vmax %d must be greater then 0.\n",
				   THISFUNCTION,options[OPT_VMAX]);
			return(-1);
		}
	}

	if (SEQUENCE) {
		*local_method = 2;
		if (*architecture == 0) {
			*ndims_tot = 1;
		}
		else if (*architecture > 0) {
			mesh_dims[0] = 2;
			mesh_dims[1] = mesh_dims[2] = 1;
		}
		*ndims = 1;
		goto End_Label;
	}

	/* Get local method, if any */
	*local_method = 0;
	if (*global_method == 1)
		*local_method = 1;
	else {
		if (PROMPT) {
			printf("Local refinement method:\n");
			printf("  (1) Kernighan-Lin\n");
			printf("  (2) None\n");
		}
		*local_method = (int)options[OPT_LOCAL];
		if (*local_method < 1 || *local_method > 2) {
			printf("%s -- Local method %d must be 1 and 2.\n",
				   THISFUNCTION,options[OPT_LOCAL]);
			return(-1);
		}
	}

	/* Now learn about the parallel architecture. */
	if (*architecture == 0) {
	/* Get total number of hypercube dimensions in which to partition. */
		*ndims_tot = 0;
		if (PROMPT)
			printf("Total number of target hypercube dimensions: ");
		*ndims_tot = nparts[0];
		if (*ndims_tot < 1) {
			printf(" Number of divisions must be at least 1\n");
			printf("%s -- Number of divisions %d must be at least 1.\n",
				   THISFUNCTION,nparts[0]);
			return(-1);
		}
		*nprocs = 1 << (*ndims_tot);
	}

	else {			/* Get dimensions of mesh. */
		mesh_dims[1] = mesh_dims[2] = 1;
		if (*architecture == 2) {
			if (PROMPT)
				printf("X and Y extent of of 2-D mesh: ");
			mesh_dims[0] = nparts[0];
			mesh_dims[1] = nparts[1];
		}
		else if (*architecture == 3) {
			if (PROMPT)
				printf("X, Y and Z extent of 3-D mesh: ");
			mesh_dims[0] = nparts[0];
			mesh_dims[1] = nparts[1];
			mesh_dims[2] = nparts[2];
		}
		else {			/* Anything else => 1-D mesh */
			if (PROMPT)
				printf("Size of 1-D mesh: ");
			mesh_dims[0] = nparts[0];
			*architecture = 1;
		}
		*nprocs = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
	}

	/* Get number of dimensions in which to partition at each level. */
	*ndims = 0;
	if (*nprocs <= 3) {
		*ndims = 1;
	}
	else if (*nprocs <= 7) {
		if (PROMPT) {
			printf("Partitioning dimension: \n");
			printf("  (1) Bisection\n");
			printf("  (2) Quadrisection\n");
		}
		*ndims = (int)options[OPT_NDIMS];
		if (*ndims < 1 || *ndims > 2) {
			printf("%s -- Ndims %d must be 1 or 2 for %d processors.\n",
				   THISFUNCTION,options[OPT_NDIMS],*nprocs);
			return(-1);
		}
	}
	else {
		if (PROMPT) {
			printf("Partitioning dimension: \n");
			printf("  (1) Bisection\n");
			printf("  (2) Quadrisection\n");
			printf("  (3) Octasection\n");
		}
		*ndims = (int)options[OPT_NDIMS];
		if (*ndims < 1 || *ndims > 3) {
			printf("%s -- Ndims %d must be between 1 and 3 for %d processors.\n",
				   THISFUNCTION,options[OPT_NDIMS],*nprocs);
			return(-1);
		}
	}
End_Label: 

	if (*global_method == 1 || *rqi_flag) {
		if (*vmax < 2 * (1 << *ndims)) {
			*vmax = 2 * (1 << *ndims);
		}
	}

	return(0);
}
