/*!\file InterpFromMeshToMesh3d.c
 * \brief: data interpolation from a list of (x,y,values) into mesh vertices
*/

#include "./InterpFromMeshToMesh3d.h"

void InterpFromMeshToMesh3dUsage(void){/*{{{*/
	_printf0_("INTERPFROMMESHTOMESH3D - interpolation from a 3d hexahedron mesh onto a list of point\n");
	_printf0_("\n");
	_printf0_("   This function is a multi-threaded mex file that interpolates a field\n");
	_printf0_("   defined on a triangular mesh onto a list of point\n");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("      data_prime=InterpFromMeshToMesh3d(index,x,y,z,data,x_prime,y_prime,z_prime,default_value);\n");
	_printf0_("\n");
	_printf0_("      index: index of the mesh where data is defined\n");
	_printf0_("      x,y,z: coordinates of the nodes where data is defined\n");
	_printf0_("      data: matrix holding the data to be interpolated onto the mesh.\n");
	_printf0_("      x_prime,y_prime,z_prime: coordinates of the points onto which we interpolate.\n");
	_printf0_("      default_value: default value if no data is found (holes).\n");
	_printf0_("      data_prime: vector of mesh interpolated data.\n");
	_printf0_("\n");
	_printf0_("   Example:\n");
	_printf0_("      load('temperature.mat');\n");
	_printf0_("      md.initialization.temperature=InterpFromMeshToMesh3d(index,x,y,z,temperature,md.mesh.x,md.mesh.y,md.mesh.z,253);\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(InterpFromMeshToMesh3d_python){

	/*input: */
	double* index_data=NULL;
	int     index_data_rows;

	double* x_data=NULL;
	double* y_data=NULL;
	double* z_data=NULL;

	int     x_data_rows;
	int     y_data_rows;
	int     z_data_rows;

	double* data=NULL; 
	int     data_rows;
	int     data_cols;

	double* x_prime=NULL;
	double* y_prime=NULL;
	double* z_prime=NULL;

	int     x_prime_rows;
	int     y_prime_rows;
	int     z_prime_rows;

	double  default_value;

	/*Intermediary*/
	int nods_data;
	int nels_data;
	int nods_prime;

	/* output: */
	IssmSeqVec<double>*  data_prime=NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	#ifdef _HAVE_MATLAB_MODULES_
	CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&InterpFromMeshToMesh3dUsage);
	#endif

	/*Input datasets: */
	FetchData(&index_data,&index_data_rows,NULL,INDEXHANDLE);
	FetchData(&x_data,&x_data_rows,NULL,XHANDLE);
	FetchData(&y_data,&y_data_rows,NULL,YHANDLE);
	FetchData(&z_data,&z_data_rows,NULL,ZHANDLE);
	FetchData(&data,&data_rows,&data_cols,DATAHANDLE);
	FetchData(&x_prime,&x_prime_rows,NULL,XPRIMEHANDLE);
	FetchData(&y_prime,&y_prime_rows,NULL,YPRIMEHANDLE);
	FetchData(&z_prime,&z_prime_rows,NULL,ZPRIMEHANDLE);
	FetchData(&default_value,DEFAULTHANDLE);

	/*some checks*/
	if (x_data_rows!=y_data_rows || x_data_rows!=z_data_rows){
		_error_("vectors x, y and z should have the same length!");
	}
	if (x_prime_rows!=y_prime_rows || x_prime_rows!=z_prime_rows){
		_error_("vectors x_prime, y_prime and z_prime should have the same length!");
	}
	/*get number of elements and number of nodes in the data*/
	nels_data=index_data_rows;
	nods_data=x_data_rows;
	nods_prime=x_prime_rows;

	/* Run core computations: */
	InterpFromMeshToMesh3dx(&data_prime,index_data,x_data,y_data,z_data,nods_data,nels_data,data,data_rows,x_prime,y_prime,z_prime,nods_prime,default_value);

	/*Write data: */
	WriteData(DATAPRIME,data_prime);

	/*end module: */
	xDelete<double>(index_data);
	xDelete<double>(x_data);
	xDelete<double>(y_data);
	xDelete<double>(z_data);
	xDelete<double>(data);
	xDelete<double>(x_prime);
	xDelete<double>(y_prime);
	xDelete<double>(z_prime);
	delete data_prime;

	/*end module: */
	MODULEEND();
}
