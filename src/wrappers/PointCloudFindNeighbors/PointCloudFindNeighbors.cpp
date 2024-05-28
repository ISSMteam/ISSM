/*! \file  PointCloudFindNeighbors
    \brief: flag points that are too near one another, within an array of point coordinates
*/

#include "./PointCloudFindNeighbors.h"

void PointCloudFindNeighborsUsage(void){/*{{{*/
	_printf_("   usage:\n");
	_printf_("   [flags]=PointCloudFindNeighbors(x,y,mindistance,multithread);\n");
	_printf_("   where:\n");
	_printf_("      x,y: list of points.\n");
	_printf_("      mindistance: minimum distance that should exist between points in the cloud.\n");
	_printf_("      multithread: run multithreaded or not. with multithreads, flags can get 1 and 2 values in duplicates.\n");
	_printf_("      flags: array of flags (flag==1 means point is within mindistance of another point)\n");
	_printf_("\n");
}/*}}}*/
WRAPPER(PointCloudFindNeighbors_python){

	int i,j;

	/* required input: */
	double *x = NULL;
	double *y = NULL;
	int     nods;
	double  mindistance;
	double  multithread;

	/* output: */
	IssmSeqVec<double> *flags = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&PointCloudFindNeighborsUsage);

	/*Fetch inputs: */
	FetchData(&x,&nods,NULL,XHANDLE);  
	FetchData(&y,NULL,NULL,YHANDLE);
	FetchData(&mindistance,MINDISTANCE);
	FetchData(&multithread,MULTITHREAD);

	/*Run core routine: */
	PointCloudFindNeighborsx(&flags,x,y,nods,mindistance,multithread);

	/* output: */
	WriteData(FLAGS,flags);

	/*end module: */
	xDelete<double>(x);
	xDelete<double>(y);
	delete flags;

	/*end module: */
	MODULEEND();
}
