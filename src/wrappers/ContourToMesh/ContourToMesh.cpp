/*! \file  ContourtoMesh
    \brief: takes an  contour file, and figures out which nodes or elements from the mesh  
    are inside this contour. 
*/

#include "./ContourToMesh.h"

void ContourToMeshUsage(void){/*{{{*/
	_printf_("CONTOURTOMESH - Flag the elements or nodes inside a contour\n");
	_printf_("\n");
	_printf_("      Usage: \n");
	_printf_("         [in_nod,in_elem]=ContourToMesh(index,x,y,contourname,interptype,edgevalue)\n");
	_printf_("\n");
	_printf_("         index,x,y: mesh triangulation.\n");
	_printf_("         contourname: name of .exp file containing the contours.\n");
	_printf_("         interptype: string definining type of interpolation ('element', or 'node').\n");
	_printf_("         edgevalue: integer (0, 1 or 2) defining the value associated to the nodes on the edges of the polygons.\n");
	_printf_("         in_nod: vector of flags (0 or 1), of size nods if interptype is set to 'node' or 'element and node', \n");
	_printf_("            or of size 0 otherwise.\n");
	_printf_("         in_elem: vector of flags (0 or 1), of size nel if interptype is set to 'element' or 'element and node', \n");
	_printf_("            or of size 0 otherwise.\n");
	_printf_("\n");
	_printf_("      Example: \n");
	_printf_("         in_nod=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','node',1)\n");
	_printf_("         in_elements=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','element',0)\n");
	_printf_("         [in_nodes,in_elements]=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','element and node',0)\n");
	_printf_("\n");
}/*}}}*/
WRAPPER(ContourToMesh_python){

	/* required input: */
	int       edgevalue;
	int       nel,nods;
	double   *index       = NULL;
	double   *x           = NULL;
	double   *y           = NULL;
	char     *interptype  = NULL;
	Contours *contours    = NULL;

	/* output: */
	double *in_nod  = NULL;
	double *in_elem = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on output arguments on the matlab side: */
	#ifdef _HAVE_MATLAB_MODULES_
	if(nlhs!=1 && nlhs!=2){
		ContourToMeshUsage();
		_error_("usage. See above");
	}
	#endif
	/*check on input arguments: */
	if(nrhs!=NRHS){
		ContourToMeshUsage();
		_error_("usage. See above");
	}

	/*Fetch inputs: */
	FetchData(&index,&nel,NULL,INDEX);
	FetchData(&x,&nods,NULL,X);
	FetchData(&y,NULL,NULL,Y);
	FetchData(&edgevalue,EDGEVALUE);
	FetchData(&contours,CONTOUR);
	FetchData(&interptype,INTERPTYPE);

	/*Run interpolation routine: */
	ContourToMeshx(&in_nod,&in_elem,index,x,y,contours,interptype,nel,nods,edgevalue);

	/* output: */
	WriteData(PLHS0,in_nod,nods);
	WriteData(PLHS1,in_elem,nel);

	/*Clean up*/
	xDelete<double>(index);
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<char>(interptype);
	delete contours;
	xDelete<double>(in_nod);
	xDelete<double>(in_elem);
	/*end module: */
	MODULEEND();
}
