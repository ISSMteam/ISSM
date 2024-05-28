/*! \file  ContourtoNodes
    \brief: takes a  contour file, and figures out which nodes  (x,y list)
*/

#include "./ContourToNodes.h"

void ContourToNodesUsage(void){/*{{{*/
	_printf_("   usage:\n");
	_printf_("   [flags]=ContourToNodes(x,y,contourname,edgevalue);\n");
	_printf_("   where:\n");
	_printf_("      x,y: list of nodes.\n");
	_printf_("      contourname: name of .exp file containing the contours, or resulting structure from call to expread.\n");
	_printf_("      edgevalue: integer (0, 1 or 2) defining the value associated to the nodes on the edges of the polygons.\n");
	_printf_("      flags: vector of flags (0 or 1), of size nods.\n");
	_printf_("\n");
}/*}}}*/
WRAPPER(ContourToNodes_python){

	/* input: */
	int       edgevalue,dim1,dim2,test1,test2;
	double   *x           = NULL;
	double   *y           = NULL;
	char     *contourname = NULL;
	Contours *contours    = NULL;

	/* output: */
	double *flags = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CHECKARGUMENTS(NLHS,NRHS,&ContourToNodesUsage);

	/*Fetch inputs: */
	FetchData(&x,&dim1,&dim2,XHANDLE);
	FetchData(&y,&test1,&test2,YHANDLE);
	FetchData(&edgevalue,EDGEVALUE);
	FetchData(&contours,CONTOUR);

	/*Some sanity checks*/
	if(dim1<1) _error_("x is empty");
	if(dim2<1) _error_("x is empty");
	if(test1!=dim1) _error_("x ans y do not have the same size");
	if(test2!=dim2) _error_("x ans y do not have the same size");

	/*Run x layer */
	ContourToNodesx(&flags,x,y,dim1*dim2,contours,edgevalue);

	/* output: */
	WriteData(FLAGS,flags,dim1,dim2);

	/*Clean up*/
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<char>(contourname);
	delete contours;
	xDelete<double>(flags);

	/*end module: */
	MODULEEND();
}
