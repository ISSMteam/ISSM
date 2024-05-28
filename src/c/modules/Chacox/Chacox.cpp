/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include "./Chacox.h"

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
	#ifdef _HAVE_CHACO_ //only works if Chaco library has been compiled in.

	extern int     Using_Main;                   /* is main routine being called?                */
	extern char   *PARAMS_FILENAME;              /* name of file with parameter updates          */
	extern double  EIGEN_TOLERANCE;              /* tolerance for eigen calculations             */
	extern int     OUTPUT_ASSIGN;                /* whether to write assignment to file          */
	extern int     DEBUG_MEMORY;                 /* debug memory allocation and freeing?         */
	extern int     DEBUG_TRACE;                  /* trace main execution path                    */
	extern int     DEBUG_PARAMS;                 /* debug flag for reading parameters            */
	extern long    RANDOM_SEED;                  /* seed for random number generators            */
	extern int     ECHO;                         /* controls amount of output                    */
	extern int     PROMPT;                       /* prompt for input or not?                     */
	extern int     PRINT_HEADERS;                /* print lines for output sections?             */
	extern int     MATCH_TYPE;                   /* matching routine to call                     */
	extern double  input_time;                   /* times data file input                        */
	extern double  start_time;                   /* time partitioning starts                     */
	FILE          *params_file;                  /* file with parameter value updates            */
	int            global_method;                /* global partitioning method                   */
	int            local_method;                 /* local partitioning method                    */
	double         eigtol;                       /* tolerance in eigenvector calculation         */
	int            ndims;                        /* dimension of recursive partitioning          */
	int            architecture;                 /* 0 => hypercube, d => d-dimensional mesh      */
	int            ndims_tot;                    /* total number of cube dimensions to divide    */
	int            mesh_dims[3];                 /* dimensions of mesh of processors             */
	long           seed;                         /* for random graph mutations                   */
	int            rqi_flag;                     /* use RQI/Symmlq eigensolver?                  */
	int            vmax;                         /* if so, how many vertices to coarsen down to? */
	char           outassignname[NAME_LENGTH];   /* assignment output file name                  */
	char           outfilename[NAME_LENGTH];     /* name of output file                          */
	char          *outassignptr;                 /* name or null pointer for output assignment   */
	char          *outfileptr;                   /* name or null pointer for output file         */
	int            nprocs;                       /* number of processors being divided into      */
	double         time;                         /* timing marker                                */
	int            flag;                         /* return code from input routines              */
	double        *smalloc();                    /* safe version of malloc                       */
	//double       seconds();                    /* returns elapsed time in seconds              */
	/*int sfree(), interface(), affirm();
	void input_queries()  , smalloc_stats(), read_params(), clear_timing();  */

	int i,tvwgt;
	double tgoal;

	if (DEBUG_TRACE > 0) {
		_printf_("<Entering main>\n");
	}

	if (PRINT_HEADERS) {
		_printf_("\n                    Chaco 2.0\n");
		_printf_("          Sandia National Laboratories\n\n");
	}

	Using_Main = TRUE;
	params_file = fopen(PARAMS_FILENAME, "r");
	if (params_file == NULL && DEBUG_PARAMS > 1) {
		printf("Parameter file `%s' not found; using default parameters.\n",PARAMS_FILENAME);
	}

	start_time = time = chaco_seconds();

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
		printf("%s -- Applying weights for %d vertices.\n",__FUNCT__,nvtxs);
		tvwgt = 0;
		for (i=0; i<nvtxs; i++)
			tvwgt += vwgts[i];
	}
	else {
		tvwgt = nvtxs;
		if      ( (int)options[OPT_VWGTS] && !vwgts)
			printf("%s -- Vertex weight flag=%d, but no vertex weights specified.\n",__FUNCT__,(int)options[OPT_VWGTS]);
		else if (!(int)options[OPT_VWGTS] &&  vwgts)
			printf("%s -- Vertex weight flag=%d, so specified vertex weights ignored.\n",__FUNCT__,(int)options[OPT_VWGTS]);
	}

	if ((int)options[OPT_EWGTS] && ewgts) {
		printf("%s -- Applying weights for %d edges.\n",
			   __FUNCT__,start[nvtxs]/2);
	}
	else {
		if      ( (int)options[OPT_EWGTS] && !ewgts)
			printf("%s -- Edge weight flag=%d, but no edge weights specified.\n",__FUNCT__,(int)options[OPT_EWGTS]);
		else if (!(int)options[OPT_EWGTS] &&  ewgts)
			printf("%s -- Edge weight flag=%d, so specified edge weights ignored.\n",__FUNCT__,(int)options[OPT_EWGTS]);
	}

    if (goal) {
        printf("%s -- Applying goals for %d sets.\n",
               __FUNCT__,nprocs);
        tgoal = 0.;
        for (i=0; i<nprocs; i++)
            tgoal += goal[i];
        for (i=0; i<nprocs; i++)
            goal[i] *= (double)tvwgt/tgoal;
    }

	input_time += chaco_seconds() - time;

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

	printf("\n%s -- Calling Chaco interface:\n\n",__FUNCT__);
	flag = interface(nvtxs, start, adjacency,
		  ((int)options[OPT_VWGTS] && vwgts ? vwgts : NULL),
		  ((int)options[OPT_EWGTS] && ewgts ? ewgts : NULL),
		  x, y, z,
		  outassignptr, outfileptr,
		  assignment,
		  architecture, ndims_tot, mesh_dims, goal,
		  global_method, local_method, rqi_flag, vmax, ndims,
		  eigtol, seed);
	printf("\n%s -- Chaco interface returning flag=%d.\n",__FUNCT__,flag);

/*  Reset adjacency matrix in case calling function needs it.  */

	for (i=0; i<start[nvtxs]; adjacency[i++]--);

	if (DEBUG_MEMORY > 0) {
		_printf_("\n");
		smalloc_stats();
	}

	if (params_file != NULL)
		fclose(params_file);

	if (DEBUG_TRACE > 1) {
		_printf_("<Leaving main>\n");
	}

	return(0);

	#else //ifdef _HAVE_CHACO_
	return (0);
	#endif
}
