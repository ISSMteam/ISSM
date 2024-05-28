/*\file DistanceToMaskBoundary.c
 *\brief: compute distance from any point in a mesh to a mask boundary
 */

#include "./DistanceToMaskBoundary.h"

void DistanceToMaskBoundaryUsage(void){/*{{{*/
	_printf0_("DISTANCETOMASKBOUNDARYUSAGE - compute distance from any point in a mesh to a mask boundary\n");
	_printf0_("\n");
	_printf0_("   This function is a multi-threaded mex file\n");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("      [distance]=DistanceToMaskBoundary(x,y,mask)\n");
	_printf0_("\n");
	_printf0_("      x,y,mask: mesh vertices with corresponding mask values. \n");
	_printf0_("      distance: distance from x,y to the mask transition between 0 and 1\n");
	_printf0_("\n");
}/*}}}*/

WRAPPER(DistanceToMaskBoundary_python){

	/*input datasets: */
	double* x=NULL;
	double* y=NULL;
	double* mask=NULL;
	int     nods;

	/* output datasets: */
	double* distance=NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	#ifdef _HAVE_MATLAB_MODULES_
	CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&DistanceToMaskBoundaryUsage);
	#endif

	/*Input datasets: */
	FetchData(&x,&nods,NULL,X);
	FetchData(&y,NULL,NULL,Y);
	FetchData(&mask,NULL,NULL,MASK);

	/*Call core of computation: */
	_error_("messing up with AD, is this fonction really used??");
	//DistanceToMaskBoundaryx(&distance,x,y,mask,nods);

	/*Write results: */
	WriteData(DISTANCE,distance,nods);

	/*Free resources: */
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<double>(mask);
	xDelete<double>(distance);

	/*end module: */
	MODULEEND();
}
