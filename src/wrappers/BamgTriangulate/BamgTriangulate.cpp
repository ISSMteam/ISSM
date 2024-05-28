/*\file BamgTriangulate.c
 *\brief: bamg module.
 */
#include "./BamgTriangulate.h"

void BamgTriangulateUsage(void){/*{{{*/
	_printf0_("BAMGTRIANGULATE - Delaunay Triangulation of a list of points");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("      index=BamgTriangulate(x,y);\n");
	_printf0_("      index: index of the triangulation\n");
	_printf0_("      x,y: coordinates of the nodes\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(BamgTriangulate_python){

	/*input: */
	double* x=NULL;
	double* y=NULL;
	int     x_cols;
	int     y_rows,y_cols;
	int     nods;

	/*Output*/
	int* index=NULL;
	int  nels;

	/*Intermediary*/
	int verbose=0;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	/* CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&BamgTriangulateUsage); */
	CHECKARGUMENTS(NLHS,NRHS,&BamgTriangulateUsage);

	/*Input datasets: */
	if (verbose) _printf_("Fetching inputs\n");
	FetchData(&x,&nods,&x_cols,XHANDLE);
	FetchData(&y,&y_rows,&y_cols,YHANDLE);

	/*Check inputs*/
	if(y_rows!=nods)         _error_("x and y do not have the same length");
	if(x_cols>1 || y_cols>1) _error_("x and y should have only one column");
	if(nods<3)               _error_("At least 3 points are required");

	/* Run core computations: */
	if (verbose) _printf_("Call core\n");
	BamgTriangulatex(&index,&nels,x,y,nods);

	/*Write output*/
	WriteData(INDEX,index,nels,3);

	/*Clean up*/
	xDelete<int>(index);
	xDelete<double>(x);
	xDelete<double>(y);

	/*end module: */
	MODULEEND();
}
