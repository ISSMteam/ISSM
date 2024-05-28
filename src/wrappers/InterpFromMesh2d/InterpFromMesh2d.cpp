/*!\file InterpFromMesh2d.c
 * \brief: data interpolation from a list of (x,y,values) into mesh vertices
*/

#include "./InterpFromMesh2d.h"

void InterpFromMesh2dUsage(void){/*{{{*/
	_printf0_("   usage:\n");
	_printf0_("         data_prime=InterpFromMesh2d(index,x,y,data,x_prime,y_prime);\n");
	_printf0_("      or data_prime=InterpFromMesh2d(index,x,y,data,x_prime,y_prime,default_value);\n");
	_printf0_("      or data_prime=InterpFromMesh2d(index,x,y,data,x_prime,y_prime,default_value,contourname);\n");
	_printf0_("   where:\n");
	_printf0_("      x,y: coordinates of the nodes where data is defined\n");
	_printf0_("      index: index of the mesh where data is defined\n");
	_printf0_("      data - vector holding the data to be interpolated onto the points.\n");
	_printf0_("      x_prime,y_prime: coordinates of the mesh vertices onto which we interpolate.\n");
	_printf0_("      default_value: a scalar or vector of size length(x_prime).\n");
	_printf0_("      contourname: linear interpolation will happen on all x_interp,y_interp inside the contour, default value will be adopted on the rest of the mesh.\n");
	_printf0_("      data_prime:  vector of prime interpolated data.\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(InterpFromMesh2d_python){

	/*input: */
	int*    index_data=NULL;
	int     index_data_rows;
	int     dummy;
	double* x_data=NULL;
	int     x_data_rows;
	double* y_data=NULL;
	int     y_data_rows;
	double* data=NULL; 
	int     data_rows;
	int     data_cols;
	double* x_prime=NULL;
	double* y_prime=NULL;
	int     x_prime_rows;
	int     y_prime_rows;
	double* default_values=NULL;
	int     num_default_values=0;

	/*contours*/
	int i;
	#ifdef _HAVE_MATLAB_MODULES_
	mxArray *matlabstructure = NULL;
	#endif
	Contour<double> **contours=NULL;
	int numcontours;
	Contour<double> *contouri=NULL;

	/*Intermediary*/
	int nods_data;
	int nels_data;
	int nods_prime;

	/* output: */
	IssmSeqVec<double> *data_prime = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	#ifdef _HAVE_MATLAB_MODULES_
	if(nlhs!=NLHS){
		InterpFromMesh2dUsage();
		_error_("InterpFromMeshToMesh2dUsage usage error");
	}
	#endif
	if((nrhs!=6) && (nrhs!=7) && (nrhs!=8)){
		InterpFromMesh2dUsage();
		_error_("InterpFromMeshToMesh2dUsage usage error");
	}

	/*Input datasets: */
	FetchData(&index_data,&index_data_rows,&dummy,INDEXHANDLE);
	FetchData(&x_data,&x_data_rows,NULL,XHANDLE);
	FetchData(&y_data,&y_data_rows,NULL,YHANDLE);
	FetchData(&data,&data_rows,&data_cols,DATAHANDLE);
	FetchData(&x_prime,&x_prime_rows,NULL,XPRIMEHANDLE);
	FetchData(&y_prime,&y_prime_rows,NULL,YPRIMEHANDLE);

	if(nrhs>=7){
		/*default values: */
		FetchData(&default_values,&num_default_values,DEFAULTHANDLE);
	}
	else{
		default_values=NULL;
		num_default_values=0;
	}

	if(nrhs==8){

		#ifdef _HAVE_MATLAB_MODULES_
		/*Call expread on filename to build a contour array in the matlab workspace: */
		mexCallMATLAB( 1, &matlabstructure, 1, (mxArray**)&FILENAME, "expread");

		/*contours: */
		numcontours=mxGetNumberOfElements(matlabstructure);
		contours=xNew<Contour<double> *>(numcontours);
		for(i=0;i<numcontours;i++){
			//allocate
			contouri=new Contour<double>();
			//retrieve dimension of this contour.
			contouri->nods=(int)mxGetScalar(mxGetField(matlabstructure,i,"nods"));
			//set pointers.
			contouri->x=mxGetPr(mxGetField(matlabstructure,i,"x"));
			contouri->y=mxGetPr(mxGetField(matlabstructure,i,"y"));
			*(contours+i)=contouri;
		}
		#else
		_error_("not supported");
		#endif
	}
	else{
		numcontours=0;
		contours=NULL;
	}

	/*some checks*/
	if (x_data_rows!=y_data_rows){
		_error_("vectors x and y should have the same length!");
	}
	if (x_prime_rows!=y_prime_rows){
		_error_("vectors x_prime and y_prime should have the same length!");
	}

	/*get number of elements and number of nodes in the data*/
	nels_data=index_data_rows;
	nods_data=x_data_rows;
	nods_prime=x_prime_rows;

	/* Run core computations: */
	InterpFromMesh2dx(&data_prime,index_data,x_data,y_data,nods_data,nels_data,data,data_rows,x_prime,y_prime,nods_prime,default_values,num_default_values,contours,numcontours);

	/*Write data: */
	WriteData(DATAPRIME,data_prime);

	/*end module: */
	xDelete<int>(index_data);
	xDelete<double>(x_data);
	xDelete<double>(y_data);
	xDelete<double>(data);
	xDelete<double>(x_prime);
	xDelete<double>(y_prime);
	xDelete<double>(default_values);
	delete data_prime;

	/*end module: */
	MODULEEND();
}
