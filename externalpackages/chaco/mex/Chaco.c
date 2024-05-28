
#include <stdio.h>
#include <string.h>    /*  strcasecmp  */
#include <time.h>      /*  clock,time,difftime  */
#include "mex.h"


#define THISFUNCTION "Chaco"

/* Input Arguments */

#define    A_IN         prhs[0]
#define    VWGTS_IN     prhs[1]
#define    EWGTS_IN     prhs[2]
#define    XYZ_IN       prhs[3]
#define    OPTNS_IN     prhs[4]
#define    NPARTS_IN    prhs[5]
#define    GOAL_IN      prhs[6]

/* Output Arguments */

#define    ASSGN_OUT    plhs[0]


void ChacoUsage( void );

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
);


void mexFunction(
	int           nlhs,           /* number of outputs */
	mxArray       *plhs[],        /* array of pointers to output arguments */
	int           nrhs,           /* number of inputs */
	const mxArray *prhs[]         /* array of pointers to input arguments */
)
{
    int       nvtxs;		/* number of vertices in graph */
    int      *start;		/* start of edge list for each vertex */
    int      *adjacency;	/* edge list data */
    int      *vwgts=NULL;	/* weights for all vertices */
    float    *ewgts=NULL;	/* weights for all edges */
    float    *x=NULL;
    float    *y=NULL;
    float    *z=NULL;		/* coordinates for inertial method */
    short    *assignment=NULL;	/* set number of each vtx (length nvtxs+1) */
    double   options[10]={1,1,0,0,1,1,50,0,.001,7654321};
							/* architecture and partitioning options */
    int      *nparts=NULL;	/* number of parts options */
    double   *goal=NULL;	/* desired set sizes */

	int    i, nedges, nterms, ncols, ierr=0;
	double *p;
	mwIndex *mwstart, *mwadjacency;

	clock_t clock0,clock1;
	time_t  time0, time1;

	/* Check for proper number of arguments */

	if      (nrhs == 0 && nlhs == 0) {
		ChacoUsage();
		return;
	}
	else if (nrhs <  1 || nlhs >  1) {
		ChacoUsage();
		mexErrMsgTxt(" ");
	}

	clock0=clock();
	time0 =time(NULL);
	mexPrintf("\nChaco Module -- %s",ctime(&time0));

	/* Assign pointers to the various parameters */
	/* and convert to the appropriate format for Chaco */

	if (!mxIsEmpty(A_IN) && mxIsNumeric(A_IN) && mxIsSparse(A_IN)) {
		nvtxs = mxGetN(A_IN);
		mwstart = mxGetJc(A_IN);
		start = mxMalloc((mxGetN(A_IN)+1)*sizeof(int));
		for (i=0; i<(mxGetN(A_IN)+1); i++)
			start[i]= (int)mwstart[i];
		nedges = start[nvtxs];
		mwadjacency = mxGetIr(A_IN);
		adjacency = mxMalloc(mxGetNzmax(A_IN)*sizeof(int));
		for (i=0; i<mxGetNzmax(A_IN); i++)
			adjacency[i]= (int)mwadjacency[i];
	}
	else {
		mexPrintf("%s -- Adjacency matrix must be numeric and sparse.\n",
				  THISFUNCTION);
		mexErrMsgTxt(" ");
	}

	if (nrhs >= 2 && !mxIsEmpty(VWGTS_IN)) {
		if (mxIsNumeric(VWGTS_IN) && !mxIsSparse(VWGTS_IN) &&
			(nterms=mxGetM(VWGTS_IN)*mxGetN(VWGTS_IN)) == nvtxs) {
			vwgts = (int *) mxCalloc(nvtxs, sizeof(int));
			p = mxGetPr(VWGTS_IN);
			for (i = 0; i < nvtxs; vwgts[i++] = (int) *p++);
		}
		else {
			mexPrintf("%s -- Vertex weight vector must be numeric, full, and length=%d.\n",
					  THISFUNCTION,nvtxs);
			mexErrMsgTxt(" ");
		}
	}

	if (nrhs >= 3 && !mxIsEmpty(EWGTS_IN)) {
		if (mxIsNumeric(EWGTS_IN) && mxIsSparse(EWGTS_IN) &&
			(nterms=mxGetNzmax(EWGTS_IN)) == nedges) {
			ewgts = (float *) mxCalloc(nedges, sizeof(float));
			p = mxGetPr(A_IN);
			for (i = 0; i < nedges; ewgts[i++] = (float) *p++);
		}
		else {
			mexPrintf("%s -- Edge weight matrix must be numeric, sparse, and nonzeroes=%d.\n",
					  THISFUNCTION,nedges);
			mexErrMsgTxt(" ");
		}
	}

	if (nrhs >= 4 && (ncols = mxGetN(XYZ_IN)) >= 1) {
		if (mxIsNumeric(XYZ_IN) && !mxIsSparse(XYZ_IN) &&
			mxGetM(XYZ_IN) == nvtxs) {
			x = (float *) mxCalloc(nvtxs, sizeof(float));
			p = mxGetPr(XYZ_IN);
			for (i = 0; i < nvtxs; x[i++] = (float) *p++);
		}
		else {
			mexPrintf("%s -- XYZ coordinate matrix must be numeric, full, and length=%d.\n",
					  THISFUNCTION,nvtxs);
			mexErrMsgTxt(" ");
		}
		if ((ncols                 ) >= 2) {
			y = (float *) mxCalloc(nvtxs, sizeof(float));
			for (i = 0; i < nvtxs; y[i++] = (float) *p++);
			if ((ncols                 ) >= 3) {
				z = (float *) mxCalloc(nvtxs, sizeof(float));
				for (i = 0; i < nvtxs; z[i++] = (float) *p++);
			}
		}
	}

	if (nrhs >= 5 && !mxIsEmpty(OPTNS_IN)) {
		if (mxIsNumeric(OPTNS_IN) && !mxIsSparse(OPTNS_IN) &&
			(nterms=mxGetM(OPTNS_IN)*mxGetN(OPTNS_IN))) {
			p = mxGetPr(OPTNS_IN);
			for (i = 0; i < (nterms<10 ? nterms : 10); options[i++] = *p++);
		}
		else {
			mexPrintf("%s -- Options vector must be numeric and full.\n",
					  THISFUNCTION);
			mexErrMsgTxt(" ");
		}
	}

	if (nrhs >= 6 && !mxIsEmpty(NPARTS_IN)) {
		if (mxIsNumeric(NPARTS_IN) && !mxIsSparse(NPARTS_IN) &&
			(nterms=mxGetM(NPARTS_IN)*mxGetN(NPARTS_IN))) {
			nparts = (int *) mxCalloc(nterms, sizeof(int));
			p = mxGetPr(NPARTS_IN);
			for (i = 0; i < nterms; nparts[i++] = (int) *p++);
		}
		else {
			mexPrintf("%s -- Parts vector must be numeric and full.\n",
					  THISFUNCTION);
			mexErrMsgTxt(" ");
		}
	}

	if (nrhs >= 7 && !mxIsEmpty(GOAL_IN)) {
		if (mxIsNumeric(GOAL_IN) && !mxIsSparse(GOAL_IN) &&
			(nterms=mxGetM(GOAL_IN)*mxGetN(GOAL_IN))) {
			goal = (double *) mxCalloc(nterms, sizeof(double));
			p = mxGetPr(GOAL_IN);
			for (i = 0; i < nterms; goal[i++] = *p++);
		}
		else {
			mexPrintf("%s -- Goal vector must be numeric and full.\n",
					  THISFUNCTION);
			mexErrMsgTxt(" ");
		}
	}

	assignment = (short *) mxCalloc(nvtxs, sizeof(short));

    /* Do the actual computations in a subroutine */

	ierr = Chacox(nvtxs, start, adjacency, vwgts, ewgts, x, y, z,
		assignment, options, nparts, goal);

    /* Create matrices for the return arguments */

	if (!ierr) {
		ASSGN_OUT = mxCreateDoubleMatrix(1,nvtxs,mxREAL);
		p = mxGetPr(ASSGN_OUT);
		for (i = 0; i < nvtxs; *p++ = (double) assignment[i++]);
	}
	else
		ASSGN_OUT = mxCreateDoubleMatrix(0,0,mxREAL);

	/* Free what we allocated */
   
	if (!assignment) mxFree((void *) assignment);
	if (!goal)       mxFree((void *) goal);
	if (!nparts)     mxFree((void *) nparts);
	if (!z)          mxFree((void *) z);
	if (!y)          mxFree((void *) y);
	if (!x)          mxFree((void *) x);
	if (!ewgts)      mxFree((void *) ewgts);
	if (!vwgts)      mxFree((void *) vwgts);
	if (!adjacency)  mxFree((void *) adjacency);
	if (!start)      mxFree((void *) start);

    clock1=clock();
    time1 =time(NULL);
    mexPrintf("Chaco Module -- %f CPU seconds; %f elapsed seconds.\n\n",
              ((double)(clock1-clock0))/CLOCKS_PER_SEC,difftime(time1,time0));

	return;
}

void ChacoUsage( void )
{
	mexPrintf("\n");
	mexPrintf("Usage: [assgn] = Chaco(A,vwgts,ewgts,xyz,options,nparts,goal);\n");
	mexPrintf("\n");

	return;
}

