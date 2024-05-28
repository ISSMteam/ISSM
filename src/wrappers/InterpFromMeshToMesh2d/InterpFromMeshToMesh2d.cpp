/*\file InterpFromMeshToMesh2d.c
 *\brief: bamg module.
 */
#include "./InterpFromMeshToMesh2d.h"

void InterpFromMeshToMesh2dUsage(void){/*{{{*/
	_printf0_("INTERFROMMESHTOMESH2D - interpolation from a 2d triangular mesh onto a list of point\n");
	_printf0_("\n");
	_printf0_("   This function is a multi-threaded mex file that interpolates a field\n");
	_printf0_("   defined on a Delaunay triangulation onto a list of point\n");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("         data_interp=InterpFromMeshToMesh2d(index,x,y,data,x_interp,y_interp);\n");
	_printf0_("      or data_interp=InterpFromMeshToMesh2d(index,x,y,data,x_interp,y_interp,OPTIONS);\n");
	_printf0_("\n");
	_printf0_("      index             : index of the mesh where data is defined (e.g. md.mesh.elements)\n");
	_printf0_("      x,y               : coordinates of the nodes where data is defined\n");
	_printf0_("      data              : matrix holding the data to be interpolated onto the mesh. (one column per field)\n");
	_printf0_("      x_interp,y_interp : coordinates of the points onto which we interpolate.\n");
	_printf0_("      data_interp       : vector of mesh interpolated data.\n");
	_printf0_("      Available options :\n");
	_printf0_("         - 'default' : default value if point is outsite of triangulation (instead of linear interpolation)\n");
	_printf0_("\n");
	_printf0_("   Example:\n");
	_printf0_("      load('temperature.mat');\n");
	_printf0_("      md.initialization.temperature=InterpFromMeshToMesh2d(index,x,y,temperature,md.mesh.x,md.mesh.y);\n");
	_printf0_("      md.initialization.temperature=InterpFromMeshToMesh2d(index,x,y,temperature,md.mesh.x,md.mesh.y,'default',253);\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(InterpFromMeshToMesh2d_python){

	/*Intermediaties*/
	int*     index              = NULL;
	double*  x_data             = NULL;
	double*  y_data             = NULL;
	double*  data               = NULL;
	int      nods_data,nels_data;
	int      M_data,N_data;
	double*  x_interp           = NULL;
	double*  y_interp           = NULL;
	int      N_interp;
	Options* options   = NULL;
	double*  data_interp = NULL;
	int      test1,test2,test;

	/*Boot module: */
	MODULEBOOT();

	/*checks on output arguments on the matlab side: */
	#ifdef _HAVE_MATLAB_MODULES_
	if(nlhs!=NLHS){
		InterpFromMeshToMesh2dUsage();
		_error_("InterpFromMeshToMesh2dUsage usage error");
	}
	#endif
	/*check on input arguments: */
	if(nrhs<NRHS){
		InterpFromMeshToMesh2dUsage();
		_error_("InterpFromMeshToMesh2dUsage usage error");
	}

	/*Fetch inputs: */
	FetchData(&index,&nels_data,&test,INDEX); if(test!=3) _error_("index should have 3 columns");
	FetchData(&x_data,&nods_data,X);          if(nods_data<3) _error_("there should be at least three points");
	FetchData(&y_data,&test,Y);               if(test!=nods_data) _error_("vectors x and y should have the same length");
	FetchData(&data,&M_data,&N_data,DATA);    if(M_data*N_data<1) _error_("data is empty");
	FetchData(&x_interp,&N_interp,XINTERP);   if(N_interp<1) _error_("no interpolation requested");
	FetchData(&y_interp,&test,YINTERP);       if(test!=N_interp) _error_("vectors x_interp and y_interp should have the same length");
	FetchData(&options,NRHS,nrhs,ARGUMENTS);

	/*Run core computations*/
	InterpFromMeshToMesh2dx(&data_interp,index,x_data,y_data,nods_data,nels_data,data,M_data,N_data,x_interp,y_interp,N_interp,options);

	/*Write data: */
	WriteData(DATAINTERP,data_interp,N_interp,N_data);

	/*end module: */
	xDelete<int>(index);
	xDelete<double>(x_data);
	xDelete<double>(y_data);
	xDelete<double>(data);
	xDelete<double>(x_interp);
	xDelete<double>(y_interp);
	xDelete<double>(data_interp);
	delete options;
	MODULEEND();
}
