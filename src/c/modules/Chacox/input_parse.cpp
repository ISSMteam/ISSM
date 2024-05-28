/*!\file:  input_parse.cpp
 * \brief  needed by Chacox.cpp
 */ 

#include "./Chacox.h"

#undef __FUNCT__ 
#define __FUNCT__  "input_parse"

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

	#ifdef _HAVE_CHACO_ //only works if Chaco library has been compiled in.

	extern int SEQUENCE;	/* sequence instead of partition graph? */
	extern int ARCHITECTURE;	/* 0=> hypercube, d=> d-dimensional mesh */
	extern int OUTPUT_ASSIGN;	/* write assignments to file? */
	extern int ECHO;		/* copy input to screen? results to file? */
	extern int DEBUG_TRACE;	/* trace main execution path */
	extern int PROMPT;		/* prompt for input? */
	extern int MATCH_TYPE;      /* max-matching routine to call */
	int       eigensolver;	/* which kind of eigensolver to use */

	if (DEBUG_TRACE > 0) {
		_printf_("<Entering input_parse>\n");
	}

	if (PROMPT) {
		_printf_("Parallel machine architecture:\n");
		_printf_("  (0) Hypercube\n");
		_printf_("  (1) One-dimensional mesh\n");
		_printf_("  (2) Two-dimensional mesh\n");
		_printf_("  (3) Three-dimensional mesh\n");
	}
	*architecture = (int)options[OPT_ARCH];
	if (*architecture < 0 || *architecture > 3) {
		printf("%s -- Architecture %d must be between 0 and 3.\n",__FUNCT__,*architecture);
		return(-1);
	}

	/* Name output assignment file. */
	if (PROMPT)
		_printf_("Assignment output file: ");
	outassignname = NULL;

	/* Name output results file. */
	if (PROMPT)
		_printf_("File name for saving run results: ");
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
			_printf_("Global partitioning method:\n");
			_printf_("  (1) Multilevel-KL\n");
			_printf_("  (2) Spectral\n");
			_printf_("  (3) Inertial\n");
			_printf_("  (4) Linear\n");
			_printf_("  (5) Random\n");
			_printf_("  (6) Scattered\n");
			_printf_("  (7) Read-from-file\n");
		}
		*global_method = (int)options[OPT_GLOBAL];
		if (*global_method < 1 || *global_method > 7) {
			printf("%s -- Global method %d must be between 1 and 7.\n",__FUNCT__,*global_method);
			return(-1);
		}
	}

	if (*global_method == 7) {	/* Name and open input assignment file. */
		if (PROMPT)
			_printf_("Assignment input file: ");
	}

	else if (*global_method == 3) {
		if (PROMPT)
			_printf_("Geometry input file name: ");
	}

	else if (*global_method == 2) {
		if (PROMPT) {
			_printf_("Eigensolver:\n");
			_printf_("  (1) Multilevel RQI/Symmlq\n");
			_printf_("  (2) Lanczos\n"); 
		}
		eigensolver = (int)options[OPT_RQI];
		if (eigensolver < 0 || eigensolver > 2) {
			printf("%s -- RQI/Symmlq flag %d must be between 0 and 2.\n",__FUNCT__,eigensolver);
			return(-1);
		}
		if (eigensolver == 1) {
			if (MATCH_TYPE == 5) {	/* geometric matching */
				if (PROMPT)
					_printf_("Geometry input file name: ");
			}
			*rqi_flag = 1;
			if (PROMPT)
				_printf_("Number of vertices to coarsen down to: ");
			*vmax = (int)options[OPT_VMAX];
			if (*vmax <= 0) {
				printf("%s -- Vmax %d must be greater then 0.\n",__FUNCT__,*vmax);
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
				_printf_("Geometry input file name: ");
		}
		if (PROMPT)
			_printf_("Number of vertices to coarsen down to: ");
		*vmax = (int)options[OPT_VMAX];
		if (*vmax <= 0) {
			printf("%s -- Vmax %d must be greater then 0.\n",__FUNCT__,*vmax);
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
			_printf_("Local refinement method:\n");
			_printf_("  (1) Kernighan-Lin\n");
			_printf_("  (2) None\n");
		}
		*local_method = (int)options[OPT_LOCAL];
		if (*local_method < 1 || *local_method > 2) {
			printf("%s -- Local method %d must be 1 and 2.\n",__FUNCT__,*local_method);
			return(-1);
		}
	}

	/* Now learn about the parallel architecture. */
	if (*architecture == 0) {
	/* Get total number of hypercube dimensions in which to partition. */
		*ndims_tot = 0;
		if (PROMPT)
			_printf_("Total number of target hypercube dimensions: ");
		*ndims_tot = nparts[0];
		if (*ndims_tot < 1) {
			_printf_(" Number of divisions must be at least 1\n");
			printf("%s -- Number of divisions %d must be at least 1.\n",
				   __FUNCT__,nparts[0]);
			return(-1);
		}
		*nprocs = 1 << (*ndims_tot);
	}

	else {			/* Get dimensions of mesh. */
		mesh_dims[1] = mesh_dims[2] = 1;
		if (*architecture == 2) {
			if (PROMPT)
				_printf_("X and Y extent of of 2-D mesh: ");
			mesh_dims[0] = nparts[0];
			mesh_dims[1] = nparts[1];
		}
		else if (*architecture == 3) {
			if (PROMPT)
				_printf_("X, Y and Z extent of 3-D mesh: ");
			mesh_dims[0] = nparts[0];
			mesh_dims[1] = nparts[1];
			mesh_dims[2] = nparts[2];
		}
		else {			/* Anything else => 1-D mesh */
			if (PROMPT)
				_printf_("Size of 1-D mesh: ");
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
			_printf_("Partitioning dimension: \n");
			_printf_("  (1) Bisection\n");
			_printf_("  (2) Quadrisection\n");
		}
		*ndims = (int)options[OPT_NDIMS];
		if (*ndims < 1 || *ndims > 2) {
			printf("%s -- Ndims %d must be 1 or 2 for %d processors.\n",__FUNCT__,*ndims,*nprocs);
			return(-1);
		}
	}
	else {
		if (PROMPT) {
			_printf_("Partitioning dimension: \n");
			_printf_("  (1) Bisection\n");
			_printf_("  (2) Quadrisection\n");
			_printf_("  (3) Octasection\n");
		}
		*ndims = (int)options[OPT_NDIMS];
		if (*ndims < 1 || *ndims > 3) {
			printf("%s -- Ndims %d must be between 1 and 3 for %d processors.\n",__FUNCT__,*ndims,*nprocs);
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

	#else //#ifdef _HAVE_CHACO_ 
	return(0);
	#endif
}
