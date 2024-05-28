/*\file FetchMatlabData.cpp:
 *\brief: general I/O interface to fetch data in matlab
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./matlabio.h"
#include <cstring>

/*Primitive data types*/
void FetchData(double** pmatrix,int* pM,int *pN,const mxArray* dataref){/*{{{*/

	double*  outmatrix=NULL;
	int      outmatrix_rows,outmatrix_cols;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if( mxIsClass(dataref,"double") ||
				mxIsClass(dataref,"single") ||
				mxIsClass(dataref,"int16") ||
				mxIsClass(dataref,"int8") ||
				mxIsClass(dataref,"uint8")){
		/*Check dataref is not pointing to NaN: */
		if ( mxIsNaN(*(mxGetPr(dataref))) && (mxGetM(dataref)==1) && (mxGetN(dataref)==1) ){
			outmatrix_rows=0;
			outmatrix_cols=0;
			outmatrix=NULL;
		}
		else{
			if(!mxIsClass(dataref,"double") && !mxIsClass(dataref,"single")){
				_printf_("Warning: converting matlab data from '" << mxGetClassName(dataref) << "' to 'double'\n");
			}
			/*Convert matlab matrix to double* matrix: */
			MatlabMatrixToDoubleMatrix(&outmatrix,&outmatrix_rows,&outmatrix_cols,dataref);
		}
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM)*pM=outmatrix_rows;
	if (pN)*pN=outmatrix_cols;

}
/*}}}*/
void FetchData(int** pmatrix,int* pM,int *pN,const mxArray* dataref){/*{{{*/

	int     i,outmatrix_rows,outmatrix_cols;
	double *doublematrix=NULL;
	int    *outmatrix=NULL;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if( mxIsClass(dataref,"double") ||
				mxIsClass(dataref,"single") ||
				mxIsClass(dataref,"int16") ||
				mxIsClass(dataref,"int8") ||
				mxIsClass(dataref,"uint8")){

		/*Check dataref is not pointing to NaN: */
		if ( mxIsNaN(*(mxGetPr(dataref))) && (mxGetM(dataref)==1) && (mxGetN(dataref)==1) ){
			outmatrix_rows=0;
			outmatrix_cols=0;
			outmatrix=NULL;
		}
		else{
			if(!mxIsClass(dataref,"double") && !mxIsClass(dataref,"single")){
				_printf_("Warning: converting matlab data from '" << mxGetClassName(dataref) << "' to 'double'\n");
			}
			/*Convert matlab matrix to double* matrix: */
			MatlabMatrixToDoubleMatrix(&doublematrix,&outmatrix_rows,&outmatrix_cols,dataref);

			/*Convert double matrix into integer matrix: */
			outmatrix=xNew<int>(outmatrix_rows*outmatrix_cols);
			for(i=0;i<outmatrix_rows*outmatrix_cols;i++)outmatrix[i]=(int)doublematrix[i];
		}
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM)*pM=outmatrix_rows;
	if (pN)*pN=outmatrix_cols;
}
/*}}}*/
void FetchData(bool** pmatrix,int* pM,int *pN,const mxArray* dataref){/*{{{*/

	int     i,outmatrix_rows,outmatrix_cols;
	double *doublematrix=NULL;
	bool   *outmatrix=NULL;

	if(mxIsEmpty(dataref) ){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if (mxIsClass(dataref,"double") ){

		/*Check dataref is not pointing to NaN: */
		if ( mxIsNaN(*(mxGetPr(dataref))) && (mxGetM(dataref)==1) && (mxGetN(dataref)==1) ){
			outmatrix_rows=0;
			outmatrix_cols=0;
			outmatrix=NULL;
		}
		else{

			/*Convert matlab matrix to double* matrix: */
			MatlabMatrixToDoubleMatrix(&doublematrix,&outmatrix_rows,&outmatrix_cols,dataref);

			/*Convert double matrix into integer matrix: */
			outmatrix=xNew<bool>(outmatrix_rows*outmatrix_cols);
			for(i=0;i<outmatrix_rows;i++)outmatrix[i]=(bool)doublematrix[i];
		}
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM)*pM=outmatrix_rows;
	if (pN)*pN=outmatrix_cols;
}
/*}}}*/
void FetchData(double** pvector,int* pM,const mxArray* dataref){/*{{{*/

	double* outvector=NULL;
	int M,N;

	/*Use Fetch matrix*/
	FetchData(&outvector,&M,&N,dataref) ;

	/*Check that it is a vector*/
	if(M*N>0 && (M!=1 && N!=1)){
		_error_("input vector of size " << M << "x" << N << " should have only one column");
	}

	/*Transpose Row vectors*/
	if(M==1 && N>1) M=N;

	/*Assign output pointers:*/
	*pvector=outvector;
	if(pM)*pM=M;
}
/*}}}*/
void FetchData(int** pvector,int* pM,const mxArray* dataref){/*{{{*/

	int    i;
	double *doublevector   = NULL;
	int    *outvector      = NULL;
	int     outvector_rows;

	if(mxIsEmpty(dataref)){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outvector_rows=0;
		outvector=NULL;
	}
	else if (mxIsClass(dataref,"double") ){

		/*Convert matlab vector to double*  vector: */
		FetchData(&doublevector,&outvector_rows,dataref);

		/*Convert double vector into integer vector: */
		outvector=xNew<int>(outvector_rows);
		for(i=0;i<outvector_rows;i++)outvector[i]=(int)doublevector[i];
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pvector=outvector;
	if (pM)*pM=outvector_rows;
}
/*}}}*/
void FetchData(bool** pvector,int* pM,const mxArray* dataref){/*{{{*/

	int    i;
	double *doublevector   = NULL;
	bool   *outvector      = NULL;
	int     outvector_rows;

	if(mxIsEmpty(dataref)){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outvector_rows=0;
		outvector=NULL;
	}
	else if (mxIsClass(dataref,"double") ){

		/*Convert matlab vector to double*  vector: */
		FetchData(&doublevector,&outvector_rows,dataref);

		/*Convert double vector into integer vector: */
		outvector=xNew<bool>(outvector_rows);
		for(i=0;i<outvector_rows;i++)outvector[i]=(bool)doublevector[i];
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pvector=outvector;
	if (pM)*pM=outvector_rows;
}
/*}}}*/
void FetchData(float** pvector,int* pM,const mxArray* dataref){/*{{{*/

	int    i;
	double *doublevector   = NULL;
	float  *outvector      = NULL;
	int     outvector_rows;

	if(mxIsEmpty(dataref)){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outvector_rows=0;
		outvector=NULL;
	}
	else if (mxIsClass(dataref,"double") ){

		/*Convert matlab vector to double*  vector: */
		FetchData(&doublevector,&outvector_rows,dataref);

		/*Convert double vector into float vector: */
		outvector=xNew<float>(outvector_rows);
		for(i=0;i<outvector_rows;i++)outvector[i]=(float)doublevector[i];
	}
	else{
		/*This is an error: we don't have the correct input!: */
		_error_("Input parameter of class " << mxGetClassName(dataref) << " not supported yet");
	}

	/*Assign output pointers:*/
	*pvector=outvector;
	if (pM)*pM=outvector_rows;
}
/*}}}*/
void FetchData(char** pstring,const mxArray* dataref){/*{{{*/

	char* outstring=NULL;

	/*Ok, the string should be coming directly from the matlab workspace: */
	if (!mxIsClass(dataref,"char")){
		_error_("input data_type is not a string!");
	}
	else{
		/*Recover the string:*/
		int stringlen;

		stringlen = mxGetM(dataref)*mxGetN(dataref)+1;
		outstring =xNew<char>(stringlen);
		mxGetString(dataref,outstring,stringlen);
	}

	/*Assign output pointers:*/
	*pstring=outstring;
}/*}}}*/
void FetchData(double* pscalar,const mxArray* dataref){/*{{{*/

	double scalar;

	if (!mxIsClass(dataref,"double")){
		_error_("input data_type is not a double!");
	}
	else{
		/*Recover the double: */
		scalar=mxGetScalar(dataref);
	}

	/*Assign output pointers:*/
	*pscalar=scalar;
}
/*}}}*/
void FetchData(int* pinteger,const mxArray* dataref){/*{{{*/

	int integer;

	if (!mxIsClass(dataref,"double")){
		_error_("input data_type is not a scalar!");
	}
	else{
		/*Recover the double: */
		integer=(int)mxGetScalar(dataref);
	}

	/*Assign output pointers:*/
	*pinteger=integer;
}
/*}}}*/
void FetchData(bool* pboolean,const mxArray* dataref){/*{{{*/

	bool* mxbool_ptr=NULL;

	if (mxIsClass(dataref,"logical")){
		if(mxGetM(dataref)!=1) _error_("input data is not of size 1x1");
		if(mxGetN(dataref)!=1) _error_("input data is not of size 1x1");
		mxbool_ptr=mxGetLogicals(dataref);
	}
	else{
		_error_("input data_type is not a bool!");
	}

	*pboolean=*mxbool_ptr;
}
/*}}}*/

/*ISSM objects*/
void FetchData(BamgGeom** pbamggeom,const mxArray* dataref){/*{{{*/

	/*Initialize output*/
	BamgGeom* bamggeom=new BamgGeom();

	/*Fetch all fields*/
	FetchData(&bamggeom->Vertices,&bamggeom->VerticesSize[0],&bamggeom->VerticesSize[1],mxGetAssignedField(dataref,0,"Vertices"));
	FetchData(&bamggeom->Edges, &bamggeom->EdgesSize[0], &bamggeom->EdgesSize[1], mxGetAssignedField(dataref,0,"Edges"));
	FetchData(&bamggeom->Corners, &bamggeom->CornersSize[0], &bamggeom->CornersSize[1], mxGetAssignedField(dataref,0,"Corners"));
	FetchData(&bamggeom->RequiredVertices,&bamggeom->RequiredVerticesSize[0],&bamggeom->RequiredVerticesSize[1],mxGetAssignedField(dataref,0,"RequiredVertices"));
	FetchData(&bamggeom->RequiredEdges, &bamggeom->RequiredEdgesSize[0], &bamggeom->RequiredEdgesSize[1], mxGetAssignedField(dataref,0,"RequiredEdges"));
	FetchData(&bamggeom->CrackedEdges,&bamggeom->CrackedEdgesSize[0],&bamggeom->CrackedEdgesSize[1],mxGetAssignedField(dataref,0,"CrackedEdges"));
	FetchData(&bamggeom->SubDomains,&bamggeom->SubDomainsSize[0],&bamggeom->SubDomainsSize[1],mxGetAssignedField(dataref,0,"SubDomains"));

	/*Assign output pointers:*/
	*pbamggeom=bamggeom;
}
/*}}}*/
void FetchData(BamgMesh** pbamgmesh,const mxArray* dataref){/*{{{*/

	/*Initialize output*/
	BamgMesh* bamgmesh=new BamgMesh();

	/*Fetch all fields*/
	FetchData(&bamgmesh->Vertices,&bamgmesh->VerticesSize[0],&bamgmesh->VerticesSize[1],mxGetAssignedField(dataref,0,"Vertices"));
	FetchData(&bamgmesh->Edges, &bamgmesh->EdgesSize[0], &bamgmesh->EdgesSize[1], mxGetAssignedField(dataref,0,"Edges"));
	FetchData(&bamgmesh->Triangles, &bamgmesh->TrianglesSize[0], &bamgmesh->TrianglesSize[1], mxGetAssignedField(dataref,0,"Triangles"));
	FetchData(&bamgmesh->CrackedEdges,&bamgmesh->CrackedEdgesSize[0],&bamgmesh->CrackedEdgesSize[1],mxGetAssignedField(dataref,0,"CrackedEdges"));
	FetchData(&bamgmesh->VerticesOnGeomEdge,&bamgmesh->VerticesOnGeomEdgeSize[0],&bamgmesh->VerticesOnGeomEdgeSize[1],mxGetAssignedField(dataref,0,"VerticesOnGeomEdge"));
	FetchData(&bamgmesh->VerticesOnGeomVertex,&bamgmesh->VerticesOnGeomVertexSize[0],&bamgmesh->VerticesOnGeomVertexSize[1],mxGetAssignedField(dataref,0,"VerticesOnGeomVertex"));
	FetchData(&bamgmesh->EdgesOnGeomEdge, &bamgmesh->EdgesOnGeomEdgeSize[0], &bamgmesh->EdgesOnGeomEdgeSize[1], mxGetAssignedField(dataref,0,"EdgesOnGeomEdge"));
	FetchData(&bamgmesh->IssmSegments,&bamgmesh->IssmSegmentsSize[0],&bamgmesh->IssmSegmentsSize[1],mxGetAssignedField(dataref,0,"IssmSegments"));

	/*Assign output pointers:*/
	*pbamgmesh=bamgmesh;
}
/*}}}*/
void FetchData(BamgOpts** pbamgopts,const mxArray* dataref){/*{{{*/

	/*Initialize output*/
	BamgOpts* bamgopts=new BamgOpts();

	/*Fetch all fields*/
	FetchData(&bamgopts->anisomax,mxGetField(dataref,0,"anisomax"));
	FetchData(&bamgopts->cutoff,mxGetField(dataref,0,"cutoff"));
	FetchData(&bamgopts->coeff,mxGetField(dataref,0,"coeff"));
	FetchData(&bamgopts->errg,mxGetField(dataref,0,"errg"));
	FetchData(&bamgopts->gradation,mxGetField(dataref,0,"gradation"));
	FetchData(&bamgopts->Hessiantype,mxGetField(dataref,0,"Hessiantype"));
	FetchData(&bamgopts->maxnbv,mxGetField(dataref,0,"maxnbv"));
	FetchData(&bamgopts->maxsubdiv,mxGetField(dataref,0,"maxsubdiv"));
	FetchData(&bamgopts->Metrictype,mxGetField(dataref,0,"Metrictype"));
	FetchData(&bamgopts->nbjacobi,mxGetField(dataref,0,"nbjacobi"));
	FetchData(&bamgopts->nbsmooth,mxGetField(dataref,0,"nbsmooth"));
	FetchData(&bamgopts->omega,mxGetField(dataref,0,"omega"));
	FetchData(&bamgopts->power,mxGetField(dataref,0,"power"));
	FetchData(&bamgopts->verbose,mxGetField(dataref,0,"verbose"));

	FetchData(&bamgopts->Crack,mxGetField(dataref,0,"Crack"));
	FetchData(&bamgopts->KeepVertices,mxGetField(dataref,0,"KeepVertices"));
	FetchData(&bamgopts->splitcorners,mxGetField(dataref,0,"splitcorners"));

	FetchData(&bamgopts->hmin,mxGetField(dataref,0,"hmin"));
	FetchData(&bamgopts->hmax,mxGetField(dataref,0,"hmax"));
	FetchData(&bamgopts->hminVertices,&bamgopts->hminVerticesSize[0],&bamgopts->hminVerticesSize[1],mxGetField(dataref,0,"hminVertices"));
	FetchData(&bamgopts->hmaxVertices,&bamgopts->hmaxVerticesSize[0],&bamgopts->hmaxVerticesSize[1],mxGetField(dataref,0,"hmaxVertices"));
	FetchData(&bamgopts->hVertices,&bamgopts->hVerticesLength,mxGetField(dataref,0,"hVertices"));
	FetchData(&bamgopts->metric,&bamgopts->metricSize[0],&bamgopts->metricSize[1],mxGetField(dataref,0,"metric"));
	FetchData(&bamgopts->field,&bamgopts->fieldSize[0],&bamgopts->fieldSize[1],mxGetField(dataref,0,"field"));
	FetchData(&bamgopts->err,&bamgopts->errSize[0],&bamgopts->errSize[1],mxGetField(dataref,0,"err"));

	/*Additional checks*/
	bamgopts->Check();

	/*Assign output pointers:*/
	*pbamgopts=bamgopts;
}
/*}}}*/
void FetchData(Options** poptions,int istart, int nrhs,const mxArray** pdataref){/*{{{*/

	char   *name   = NULL;
	Option *option = NULL;

	/*Initialize output*/
	Options* options=new Options();

	/*Fetch all options*/
	for (int i=istart; i<nrhs; i=i+2){
		if (!mxIsClass(pdataref[i],"char")) _error_("Argument " << i+1 << " must be name of option");

		FetchData(&name,pdataref[i]);
		if(i+1 == nrhs) _error_("Argument " << i+2 << " must exist and be value of option \"" << name << "\".");

		option=(Option*)OptionParse(name,&pdataref[i+1]);
		options->AddOption(option);
		option=NULL;
	}

	/*Assign output pointers:*/
	*poptions=options;
}
/*}}}*/
void FetchData(Contours** pcontours,const mxArray* dataref){/*{{{*/

	int             numcontours,index,test1,test2;
	char            *contourname = NULL;
	Contours        *contours    = NULL;
	Contour<double> *contouri    = NULL;

	if(mxIsClass(dataref,"char")){
		FetchData(&contourname,dataref);
		contours=ExpRead<double>(contourname);
	}
	else if(mxIsClass(dataref,"struct")){

		contours=new Contours();
		numcontours=mxGetNumberOfElements(dataref);

		for(int i=0;i<numcontours;i++){

			contouri=new Contour<double>();

			index = mxGetFieldNumber(dataref,"nods");
			if(index==-1) _error_("input structure does not have a 'nods' field");
			FetchData(&contouri->nods,mxGetFieldByNumber(dataref,i,index));

			index = mxGetFieldNumber(dataref,"x");
			if(index==-1) _error_("input structure does not have a 'x' field");
			FetchData(&contouri->x,&test1,&test2,mxGetFieldByNumber(dataref,i,index));
			if(test1!=contouri->nods || test2!=1) _error_("field x should be of size ["<<contouri->nods<<" 1]");

			index = mxGetFieldNumber(dataref,"y");
			if(index==-1) _error_("input structure does not have a 'y' field");
			FetchData(&contouri->y,&test1,&test2,mxGetFieldByNumber(dataref,i,index));
			if(test1!=contouri->nods || test2!=1) _error_("field y should be of size ["<<contouri->nods<<" 1]");

			contours->AddObject(contouri);
		}
	}
	else{
		_error_("Contour is neither a string nor a structure and cannot be loaded ("<<mxGetClassName(dataref)<<" not supported)");
	}

	/*clean-up and assign output pointer*/
	xDelete<char>(contourname);
	*pcontours=contours;
}
/*}}}*/
void FetchChacoData(int* pnvtxs,int** padjacency,int** pstart,float** pewgts,const mxArray* A_IN, const mxArray* EWGTS_IN){/*{{{*/

	/*Fetch adjacency matrix: */
	int      nvtxs       = mxGetN(A_IN);
	mwIndex* mwstart     = mxGetJc(A_IN);
	mwIndex* mwadjacency = mxGetIr(A_IN);
	int      nzmax       = mxGetNzmax(A_IN);

	int* start = xNew<int>(nvtxs+1);
	for(int i=0;i<nvtxs+1;i++) start[i]=(int)mwstart[i];

	int* adjacency = xNew<int>(nzmax);
	for(int i=0; i<nzmax; i++) adjacency[i]= (int)mwadjacency[i];

	/*Get edges weights*/
	int nedges = start[nvtxs];
	float* ewgts = NULL;
	if(!mxIsEmpty(EWGTS_IN)){
		ewgts = xNewZeroInit<float>(nedges);
		double* doublepointer = mxGetPr(A_IN);
		for(int i = 0; i<nedges;i++) ewgts[i] = (float)doublepointer[i];
	}

	/*Assign output pointers*/
	*pnvtxs     = nvtxs;
	*padjacency = adjacency;
	*pstart     = start;
	*pewgts     = ewgts;
}
/*}}}*/

/*Toolkit*/
int MatlabMatrixToDoubleMatrix(double** pmatrix,int* pmatrix_rows,int* pmatrix_cols,const mxArray* mxmatrix){/*{{{*/

	/*Get Matrix size*/
	int rows=mxGetM(mxmatrix);
	int cols=mxGetN(mxmatrix);

	/*Return of Matrix is empty*/
	if(rows*cols == 0){
		*pmatrix      = NULL;
		*pmatrix_rows = rows;
		*pmatrix_cols = cols;
		return 1;
	}

   /*Initialize output*/
   double* matrix=xNewZeroInit<double>(rows*cols);

	/*First check if we are dealing with a sparse matrix: */
	if(mxIsSparse(mxmatrix)){

		/*Dealing with sparse matrix: recover size first: */
		double* pmxmatrix=(double*)mxGetPr(mxmatrix);

      /*Now, get ir,jc and pr: */
      mwIndex* ir=mxGetIr(mxmatrix);
      mwIndex* jc=mxGetJc(mxmatrix);

      /*Now, start inserting data into double* matrix: */
      int count=0;
      for(int i=0;i<cols;i++){
         for(int j=0;j<(jc[i+1]-jc[i]);j++){
            matrix[cols*ir[count]+i]=pmxmatrix[count];
            count++;
         }
      }
   }
	else if(mxIsClass(mxmatrix,"double")){
		double* pmxmatrix=(double*)mxGetPr(mxmatrix);
      for(int i=0;i<rows;i++) for(int j=0;j<cols;j++) matrix[cols*i+j]=(double)pmxmatrix[rows*j+i];
	}
	else if(mxIsClass(mxmatrix,"single")){
		float *pmxmatrix=(float*)mxGetPr(mxmatrix);
      for(int i=0;i<rows;i++) for(int j=0;j<cols;j++) matrix[cols*i+j]=(double)pmxmatrix[rows*j+i];
	}
	else if(mxIsClass(mxmatrix,"int16")){
		short int *pmxmatrix=(short*)mxGetPr(mxmatrix);
      for(int i=0;i<rows;i++) for(int j=0;j<cols;j++) matrix[cols*i+j]=(double)pmxmatrix[rows*j+i];
	}
	else if(mxIsClass(mxmatrix,"uint8")){
      char *pmxmatrix=(char*)mxGetPr(mxmatrix);
      for(int i=0;i<rows;i++) for(int j=0;j<cols;j++) matrix[cols*i+j]=(double)pmxmatrix[rows*j+i];
	}
	else{
		_error_("Matlab matrix type "<<mxGetClassName(mxmatrix)<<" Not implemented yet");
	}

	/*Assign output pointer: */
	*pmatrix=matrix;
	*pmatrix_rows=rows;
	*pmatrix_cols=cols;

	return 1;
}/*}}}*/
mxArray* mxGetAssignedField(const mxArray* pmxa_array,int number,const char* field){/*{{{*/

	/*Output*/
	mxArray *mxfield = NULL;

	if(mxIsStruct(pmxa_array)){
		mxfield = mxGetField(pmxa_array,number,field);
	}
	else{
		/*This is an object, mxGetField returns NULL in old version of matlab (we do not have access to them)*/

		/*Intermediaries*/
		mxArray    *inputs[2];
		mwSize      ndim        = 2;
		mwSize      onebyone[2] = {1,1};

		/*create index structure used in the assignment (index.type='.' and index.subs='x' for field x*/
		const char *fnames[2];
		fnames[0] = "type"; fnames[1] = "subs";
		mxArray* pindex=mxCreateStructArray( ndim,onebyone,2,fnames);
		mxSetField( pindex, 0, "type",mxCreateString("."));
		mxSetField( pindex, 0, "subs",mxCreateString(field));
		inputs[0]=(mxArray*)pmxa_array; //this is the model
		inputs[1]=pindex;

		mexCallMATLAB( 1, &mxfield, 2, (mxArray**)inputs, "subsref");
	}

	if(mxfield == NULL) _error_("Could not find field "<< field <<" in structure");

	return mxfield;
}/*}}}*/

/*Options*/
Option* OptionParse(char* name, const mxArray* prhs[]){ /*{{{*/

	/*Initialize output*/
	Option  *option = NULL;

	/*parse the value according to the matlab data type  */
	if (mxIsClass(prhs[0],"double")  && (mxGetNumberOfElements(prhs[0])==1)){
		option=(Option*)OptionDoubleParse(name,prhs);
	}
	else if(mxIsClass(prhs[0],"double")  && (mxGetNumberOfElements(prhs[0])>1)){
		option=(Option*)OptionDoubleArrayParse(name,prhs);
	}
	else if(mxIsClass(prhs[0],"char")){
		option=(Option*)OptionCharParse(name,prhs);
	}
	else {
		_error_("Second argument value of option \""<< name <<"\" is of unrecognized class \""<< mxGetClassName(prhs[0]) <<"\".");
	}

	return option;
}/*}}}*/
GenericOption<double>*  OptionDoubleParse(char* name, const mxArray* prhs[]){ /*{{{*/

	/*Initialize option*/
	GenericOption<double>* option=new GenericOption<double>();

	/*Copy name*/
	option->name =xNew<char>(strlen(name)+1);
	memcpy(option->name,name,(strlen(name)+1)*sizeof(char));

	/*Fetch data and initialize size*/
	FetchData(&option->value,prhs[0]);
   option->size[0] = 1;
	option->size[1] = 1;

	return option;
}/*}}}*/
GenericOption<double*>* OptionDoubleArrayParse(char* name, const mxArray* prhs[]){ /*{{{*/

	/*Initialize option*/
	GenericOption<double*> * option=new GenericOption<double*>();

	/*Copy name*/
	option->name =xNew<char>(strlen(name)+1);
	memcpy(option->name,name,(strlen(name)+1)*sizeof(char));

	/*check and parse the value  */
	if (!mxIsClass(prhs[0],"double")){
		_error_("Value of option \"" << option->name  << "\" must be class \"double\", not class \"" << mxGetClassName(prhs[0]) <<"\".");
	}

	/*Fetch data and initialize size*/
   FetchData(&option->value,&option->size[0],&option->size[1],prhs[0]);

	return option;
}/*}}}*/
GenericOption<char*>*   OptionCharParse(char* name, const mxArray* prhs[]){ /*{{{*/

	/*Initialize option*/
	GenericOption<char*>* option=new GenericOption<char*>();

	/*Copy name*/
	option->name =xNew<char>(strlen(name)+1);
	memcpy(option->name,name,(strlen(name)+1)*sizeof(char));

	/*check and parse the value  */
	if(!mxIsClass(prhs[0],"char")){
		_error_("Value of option \"" << option->name  << "\" must be class \"char\", not class \"" << mxGetClassName(prhs[0]) <<"\".");
	}

	/*Fetch data and initialize size*/
	FetchData(&option->value,prhs[0]);
	option->size[0] = strlen(name);
	option->size[1] = 1;

	return(option);
}/*}}}*/
