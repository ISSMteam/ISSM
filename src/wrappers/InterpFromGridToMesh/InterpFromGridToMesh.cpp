/*!\file InterpFromGridToMesh.c
 * \brief: data interpolation from a list of (x,y,values) into mesh vertices
*/

#include "./InterpFromGridToMesh.h"

void InterpFromGridToMeshUsage(void){/*{{{*/
	_printf0_("INTERPFROMGRIDTOMESH - interpolation from a grid onto a list of points\n");
	_printf0_("\n");
	_printf0_("   This function is a multi-threaded mex file that interpolates a field\n");
	_printf0_("   defined on a grid onto a list of points based on a bilinear interpolation\n");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("      data_mesh=InterpFromGridToMesh(x,y,data,x_mesh,y_mesh,default_value);\n");
	_printf0_("\n");
	_printf0_("      data: matrix holding the data to be interpolated onto the mesh.\n");
	_printf0_("      x,y: coordinates of matrix data. (x and y must be in increasing order)\n");
	_printf0_("      x_mesh,y_mesh: coordinates of the points onto which we interpolate.\n");
	_printf0_("      default_value: default value if no data is found (holes).\n");
	_printf0_("      data_mesh: vector of mesh interpolated data.\n");
	_printf0_("\n");
	_printf0_("   Example:\n");
	_printf0_("      load('velocities.mat');\n");
	_printf0_("      md.inversion.vx_obs=InterpFromGridToMesh(x_n,y_m,vx,md.mesh.x,md.mesh.y,0);\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(InterpFromGridToMesh_python){

	int i,j;

	/*input: */
	double *x = NULL;
	double *y = NULL;
	int     x_rows,y_rows;
	double *data  = NULL;
	int     data_rows,data_cols;
	double *x_mesh = NULL;
	double *y_mesh = NULL;
	int     x_mesh_rows,y_mesh_rows;
	double  default_value;
	char*   interpolationtype = NULL;

	/* output: */
	IssmSeqVec<double>*  data_mesh=NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	if(nrhs!=6 && nrhs!=7){
		InterpFromGridToMeshUsage();
		_error_("Wrong usage. See above");
	}

	/*Input datasets: */
	FetchData(&x,&x_rows,XHANDLE);
	FetchData(&y,&y_rows,YHANDLE);
	FetchData(&data,&data_rows,&data_cols,DATAHANDLE);
	FetchData(&x_mesh,&x_mesh_rows,XMESHHANDLE);
	FetchData(&y_mesh,&y_mesh_rows,YMESHHANDLE);
	FetchData(&default_value,DEFAULTHANDLE);

	/* Run core computations: */
	if(nrhs==7){
		FetchData(&interpolationtype,INTERPENUM);
		InterpFromGridToMeshx(&data_mesh, x, x_rows,  y, y_rows, data, data_rows,data_cols, x_mesh, y_mesh, x_mesh_rows,default_value,interpolationtype);
		xDelete<char>(interpolationtype);
	}
	else{
		InterpFromGridToMeshx(&data_mesh, x, x_rows,  y, y_rows, data, data_rows,data_cols, x_mesh, y_mesh, x_mesh_rows,default_value,"bilinear");
	}

	/*Write data: */
	WriteData(DATAMESH,data_mesh);

	/*end module: */
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<double>(data);
	xDelete<double>(x_mesh);
	xDelete<double>(y_mesh);
	delete data_mesh;
	MODULEEND();
}
