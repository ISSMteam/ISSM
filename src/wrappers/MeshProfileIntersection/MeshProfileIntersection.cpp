/*! \file  MeshProfileIntersection.cpp
    \brief: takes a  .exp file (made of several profiles), and figures out its intersection 
	with a mesh.

	usage:
	[segments]=MeshProfileIntersection(index,x,y,filename);

	where:
	input:
		index,x,y is a triangulation
		filename: name of Argus style .exp file containing the segments (can be groups of disconnected segments)
	output:
		segments: array made of x1,y1,x2,y2,element_id lines (x1,y1) and (x2,y2) are segment extremities for a segment 
		belonging to the elemnt_id element. there are as many lines in segments as there are segments intersecting the 
		mesh.
*/

#include "./MeshProfileIntersection.h"

void MeshProfileIntersectionUsage(void){/*{{{*/
	_printf_("   usage:\n");
	_printf_("   [segments]=MeshProfileIntersection(index,x,y,filename);\n");
	_printf_("   where:\n");
	_printf_("   input:\n");
	_printf_("        index,x,y is a triangulation\n");
	_printf_("        filename: name of Argus style .exp file containing the segments (can be groups of disconnected segments)\n");
	_printf_("   output:\n");
	_printf_("        segments: array made of x1,y1,x2,y2,element_id lines (x1,y1) and (x2,y2) are segment extremities for a segment \n");
	_printf_("        belonging to the elemnt_id element. there are as many lines in segments as there are segments intersecting the \n");
	_printf_("        mesh.\n");
}/*}}}*/
WRAPPER(MeshProfileIntersection_python){

	int i,j;

	/* required input: */
	//mesh
	double *double_index = NULL;
	int    *index        = NULL;
	int     nel;
	double *x            = NULL;
	double *y            = NULL;
	int     nods;
	int     dummy;

	//contours
	Contours         *domain      = NULL;
	Contour<double> **contours=NULL;
	int               numcontours;
	Contour<double>  *contouri=NULL;

	/* output: */
	double* segments=NULL;
	int     numsegs;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments: */
	CHECKARGUMENTS(NLHS,NRHS,&MeshProfileIntersectionUsage);

	/*Fetch inputs: */
	//index
	FetchData(&double_index,&nel,&dummy,INDEX);
	if(dummy!=3 && dummy!=6)_error_("element triangulation should be of 3 or 6 column width!");
	index=xNew<int>(nel*3);
	for(i=0;i<nel;i++){
		for(j=0;j<3;j++){
			*(index+3*i+j)=(int)*(double_index+dummy*i+j)-1; //"C" style indexing
		}
	}
	//x and y
	FetchData(&x,&nods,X);
	FetchData(&y,&dummy,Y);

	//contours
	FetchData(&domain,FILENAME);
	// MeshProfileIntersectionx should be modified to take DataSet directly (and perhaps IssmDenseMat and IssmSeqVec).
	numcontours=domain->Size();
	contours=xNew<Contour<double>*>(numcontours);
	for(i=0;i<numcontours;i++)
		*(contours+i)=(Contour<double>*)domain->GetObjectByOffset(i);

	/*Run interpolation routine: */
	MeshProfileIntersectionx(&segments,&numsegs,index,x,y,nel,nods,contours,numcontours);

	/* output: */
	WriteData(SEGMENTS,segments,numsegs,5);

	/*end module: */
	xDelete<double>(double_index);
	xDelete<int>(index);
	xDelete<double>(x);
	xDelete<double>(y);
	delete domain;
	delete contouri;
	xDelete<double>(segments);
	MODULEEND();

}
