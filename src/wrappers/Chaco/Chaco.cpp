/*\file Chaco.c
 *\brief:  Chaco partitioner mex module
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Chaco.h"

void ChacoUsage(void){/*{{{*/
	_printf0_("\n");
	_printf0_("Usage: [assgn] = Chaco(A,vwgts,ewgts,x,y,z,options,nparts,goal);\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(Chaco_python){

	int i;
	int nterms;

	/*Inputs: */
	int     nvtxs;               /* number of vertices in graph           */
	int    *start = NULL;        /* start of edge list for each vertex    */
	int    *adjacency = NULL;    /* edge list data                        */
	int    *vwgts       = NULL;  /* weights for all vertices              */
	float  *ewgts       = NULL;  /* weights for all edges                 */
	float  *x           = NULL;
	float  *y           = NULL;
	float  *z           = NULL;  /* coordinates for inertial method       */
	double  options[10] = {1,1,0,0,1,1,50,0,0.001,7654321}; /* architecture and partitioning options */
	double *in_options  = NULL;
	int     npart;
	double *goal        = NULL;   /* desired set sizes                     */

	#ifndef _HAVE_CHACO_ //only works if dakota library has been compiled in.
	_error_("Chaco not available! Cannot carry out Chaco partitioning!");
	#endif

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CHECKARGUMENTS(NLHS,NRHS,&ChacoUsage);

	/*Fetch Data*/
	FetchChacoData(&nvtxs,&adjacency,&start,&ewgts,A_IN,EWGTS_IN);
	FetchData(&vwgts,&nterms,VWGTS_IN);
	FetchData(&x,&nterms,X_IN);
	FetchData(&y,&nterms,Y_IN);
	FetchData(&z,&nterms,Z_IN);
	FetchData(&in_options,&nterms,OPTNS_IN);
	for (i=0;i<(nterms<10?nterms:10);i++) options[i]=in_options[i]; //copy in_options into default options
	FetchData(&npart,NPARTS_IN);
	//int * nparts=xNew<int>(1); nparts[0]=npart; //weird Chacox interface ain't it?
	FetchData(&goal,&nterms,GOAL_IN);

	/*Allocate output: */
	short* assignment = xNewZeroInit<short>(nvtxs);

    /*Call core: */
	Chacox(nvtxs, start, adjacency, vwgts, ewgts, x, y, z, assignment, options,&npart, goal);

    /*Output data: */
	WriteData(ASSGN_OUT,assignment,nvtxs);

	/*Free resources:*/
	xDelete<short>(assignment); 
	xDelete<double>(goal);
	//xDelete<int>(nparts);
	xDelete<float>(z);
	xDelete<float>(y);
	xDelete<float>(x);
	xDelete<float>(ewgts);
	xDelete<int>(vwgts);
	xDelete<int>(adjacency);
	xDelete<int>(start);

	/*end module: */
	MODULEEND();
}
