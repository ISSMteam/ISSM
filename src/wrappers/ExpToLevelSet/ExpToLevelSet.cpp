/*! \file  ContourtoMesh
    \brief: takes a  contour file, a cloud of points, and figures out a levelset dependent on the distance between the contour and 
	the cloud.
*/

#include "./ExpToLevelSet.h"

void ExpToLevelSetUsage(void){/*{{{*/
	_printf_("EXPTOLEVELSET - determien levelset distance between a contour and a cloud of points\n");
	_printf_("\n");
	_printf_("      Usage: \n");
	_printf_("         distance=ExpToLevelSet(x,y,contourname)\n");
	_printf_("\n");
	_printf_("         x,y: cloud point.\n");
	_printf_("         contourname: name of .exp file containing the contours.\n");
	_printf_("         distance: distance vector representing a levelset where the 0 level is one the contour segments', \n");
	_printf_("\n");
	_printf_("      Example: \n");
	_printf_("         distance=ExpToLevelSet(md.mesh.x,md.mesh.y,'Contour.exp')\n");
	_printf_("\n");
}/*}}}*/
WRAPPER(ExpToLevelSet_python){

	/*diverse: */
	int i;

	/* required input: */
	int       nods;
	double   *x          = NULL;
	double   *y          = NULL;
	char     *interptype = NULL;
	double   *flags      = NULL;
	Contours *contours   = NULL;

	/* output: */
	double *distance  = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*check on input arguments: */
	if(nrhs!=NRHS){
		ExpToLevelSetUsage();
		_error_("usage. See above");
	}

	/*Fetch inputs: */
	FetchData(&x,&nods,NULL,X);
	FetchData(&y,NULL,NULL,Y);
	FetchData(&contours,CONTOUR);

	/*Run interpolation routine: */
	ExpToLevelSetx( &distance,x,y,nods,contours);
	ContourToNodesx(&flags,x,y,nods,contours,2);

	/*Make flags into a sign, left or right, or nill: */
	for(i=0;i<nods;i++){
		if (flags[i]==0) flags[i]=-1;
		else if (flags[i]==2) flags[i]=0;
	}

	/*Multiply flags and distance: */
	for(i=0;i<nods;i++)distance[i]*=flags[i];

	/* output: */
	WriteData(PLHS0,distance,nods);

	/*Clean up*/
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<char>(interptype);
	delete contours;
	xDelete<double>(distance);
	xDelete<double>(flags);
	
	/*end module: */
	MODULEEND();
}
