/*
 * Triangle: mesh a domain using an .exp file
 */

#include "./Triangle.h"

void TriangleUsage(void){/*{{{*/
	_printf_("\n");
	_printf_("   usage: [index,x,y,segments,segmentmarkers]=Triangle(domainoutlinefilename,rifts,area) \n");
	_printf_("      where: index,x,y defines a triangulation, segments is an array made \n");
	_printf_("      of exterior segments to the mesh domain outline, segmentmarkers is an array flagging each segment, \n");
	_printf_("      outlinefilename an Argus domain outline file, \n");
	_printf_("      area is the maximum area desired for any element of the resulting mesh, \n");
	_printf_("\n");
}/*}}}*/
WRAPPER(Triangle_python){
	
	/*intermediary: */
	double    area;
	Contours *domain = NULL;
	Contours *rifts  = NULL;

	/* output: */
	int    *index             = NULL;
	double *x                 = NULL;
	double *y                 = NULL;
	int    *segments          = NULL;
	int    *segmentmarkerlist = NULL;
	int     nel,nods,nsegs;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments: */
	CHECKARGUMENTS(NLHS,NRHS,&TriangleUsage);

	/*Fetch data needed for meshing: */
	FetchData(&domain,DOMAINOUTLINE);
	FetchData(&rifts,RIFTSOUTLINE);
	FetchData(&area,AREA);

	/*call x core: */
	Trianglex(&index,&x,&y,&segments,&segmentmarkerlist,&nel,&nods,&nsegs,domain,rifts,area);

	/*write outputs: */
	WriteData(INDEX,index,nel,3);
	WriteData(X,x,nods);
	WriteData(Y,y,nods);
	WriteData(SEGMENTS,segments,nsegs,3);
	WriteData(SEGMENTMARKERLIST,segmentmarkerlist,nsegs);

	/*free resources: */
	delete domain;
	delete rifts;
	xDelete<int>(index);
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<int>(segments);
	xDelete<int>(segmentmarkerlist);

	/*end module: */
	MODULEEND();
}
