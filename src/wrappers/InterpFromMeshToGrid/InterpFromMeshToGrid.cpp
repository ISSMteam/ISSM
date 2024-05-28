/*\file InterpFromMeshToGrid.c
 *\brief: compute diff between observed and modeled velocity
 */

#include "./InterpFromMeshToGrid.h"

void InterpFromMeshToGridUsage(void){/*{{{*/
	_printf0_("INTERPFROMMESHTOGRID - interpolation of a data defined on a mesh onto a grid\n");
	_printf0_("\n");
	_printf0_("   This function is a multi-threaded mex file that interpolates a field\n");
	_printf0_("   defined on a triangular mesh onto a regular grid\n");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("      grid=InterpFromMeshToGrid(index,x,y,data,x_grid,y_grid,default_value)\n");
	_printf0_("\n");
	_printf0_("      index,x,y: delaunay triangulation defining the mesh.\n");
	_printf0_("      meshdata: vertex values of data to be interpolated.\n");
	_printf0_("      xgrid,ygrid: parameters that define the grid\n");
	_printf0_("      default_value: value of points located out of the mesh.\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(InterpFromMeshToGrid_python){

	/*inputs */
	int*    index=NULL;
	double* x=NULL;
	double* y=NULL;
	int     nel,nods;
	double* meshdata=NULL;
	int     meshdata_length;
	double* xgrid=NULL;
	double* ygrid=NULL;
	int     nlines,ncols,test;
	double  default_value;

	/* outputs */
	double* griddata=NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	#ifdef _HAVE_MATLAB_MODULES_
	CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&InterpFromMeshToGridUsage);
	#endif

	/*Input datasets: */
	FetchData(&index,&nel,&test,INDEX);
	if(test!=3) _error_("size not supported yet");
	FetchData(&x,&nods,X);
	FetchData(&y,&test,Y);
	if(test!=nods) _error_("size not supported yet");
	FetchData(&meshdata,&meshdata_length,MESHDATA);
	FetchData(&xgrid,&ncols,XGRID); 
	FetchData(&ygrid,&nlines,YGRID);
	FetchData(&default_value,DEFAULTVALUE);

	/*Call core of computation: */
	InterpFromMeshToGridx(&griddata,index,x,y,nods,nel,meshdata,meshdata_length,xgrid,ygrid,nlines,ncols,default_value);

	/*Write results: */
	WriteData(GRIDDATA,griddata,nlines,ncols);

	/*Free resources: */
	xDelete<int>(index);
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<double>(meshdata);
	xDelete<double>(griddata);
	xDelete<double>(xgrid);
	xDelete<double>(ygrid);

	/*end module: */
	MODULEEND();
}
