/*\file WriteMatlabData.c:
 *\brief: general I/O interface to write data in matlab
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./matlabio.h"
#include "./../../../c/datastructures/datastructures.h"

/*Primitive data types*/
/*FUNCTION WriteData(mxArray** pdataref,double* matrix, int M,int N){{{*/
void WriteData(mxArray** pdataref,double* matrix, int M,int N){

	mxArray *dataref  = NULL;
	double  *tmatrix  = NULL;

	if(matrix){
		/*create the matlab matrixwith Matlab's memory manager */   
		tmatrix=(double*)mxMalloc(M*N*sizeof(double));
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				tmatrix[j*M+i]=matrix[i*N+j];
			}
		}
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,(mwSize)M);
		mxSetN(dataref,(mwSize)N);
		mxSetPr(dataref,(double*)tmatrix);
	}
	else{
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
	}
	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,int* matrix, int M,int N){{{*/
void WriteData(mxArray** pdataref,int* matrix, int M,int N){

	mxArray* dataref = NULL;
	double*  tmatrix = NULL;

	if(matrix){
		/*convert to double matrix using Matlab's memory manager*/
		double* tmatrix=(double*)mxMalloc(M*N*sizeof(double));
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				tmatrix[j*M+i]=(double)matrix[i*N+j];
			}
		}
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,(mwSize)M);
		mxSetN(dataref,(mwSize)N);
		mxSetPr(dataref,(double*)tmatrix);
	}
	else{
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
	}
	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,double* vector, int M){{{*/
void WriteData(mxArray** pdataref,double* vector, int M){

	mxArray* dataref       = NULL;
	double*  vector_matlab = NULL;

	if(vector){

		/*create the matlab vector with Matlab's memory manager */
		vector_matlab=(double*)mxMalloc(M*sizeof(double));
		for(int i=0;i<M;i++) vector_matlab[i]=vector[i];
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,(mwSize)M);
		mxSetN(dataref,(mwSize)1);
		mxSetPr(dataref,vector_matlab);
	}
	else{
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
	}

	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,int* vector, int M){{{*/
void WriteData(mxArray** pdataref,int* vector, int M){

	mxArray* dataref       = NULL;
	double*  vector_matlab = NULL;

	if(vector){

		/*create the matlab vector with Matlab's memory manager */
		vector_matlab=(double*)mxMalloc(M*sizeof(double));
		for(int i=0;i<M;i++) vector_matlab[i]=double(vector[i]);
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,(mwSize)M);
		mxSetN(dataref,(mwSize)1);
		mxSetPr(dataref,vector_matlab);
	}
	else{
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
	}

	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,short* vector, int M){{{*/
void WriteData(mxArray** pdataref,short* vector, int M){

	mxArray* dataref       = NULL;
	double*  vector_matlab = NULL;

	if(vector){

		/*create the matlab vector with Matlab's memory manager */
		vector_matlab=(double*)mxMalloc(M*sizeof(double));
		for(int i=0;i<M;i++) vector_matlab[i]=double(vector[i]);
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,(mwSize)M);
		mxSetN(dataref,(mwSize)1);
		mxSetPr(dataref,vector_matlab);
	}
	else{
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
	}

	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,double scalar){{{*/
void WriteData(mxArray** pdataref,double scalar){

	*pdataref=mxCreateDoubleScalar(scalar);
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,int integer){{{*/
void WriteData(mxArray** pdataref,int integer){

		*pdataref=mxCreateDoubleScalar((double)integer);

}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,int boolean){{{*/
void WriteData(mxArray** pdataref,bool boolean){

	*pdataref=mxCreateDoubleScalar((double)boolean);

}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,char* string){{{*/
void WriteData(mxArray** pdataref,const char* string){

		*pdataref=mxCreateString(string);
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref){{{*/
void WriteData(mxArray** pdataref){

		;

}
/*}}}*/

/*ISSM objects*/
/*FUNCTION WriteData(mxArray** pdataref,BamgGeom* bamggeom){{{*/
void WriteData(mxArray** pdataref,BamgGeom* bamggeom){

	/*Intermediary*/
	int         i;
	mxArray    *dataref           = NULL;
	const int   numfields         = 8;
	const char *fnames[numfields];
	mwSize      ndim              = 2;
	mwSize      dimensions[2]     = {1,1};

	/*Initialize field names*/
	i=0;
	fnames[i++] = "Vertices";
	fnames[i++] = "Edges";
	fnames[i++] = "TangentAtEdges";
	fnames[i++] = "Corners";
	fnames[i++] = "RequiredVertices";
	fnames[i++] = "RequiredEdges";
	fnames[i++] = "CrackedEdges";
	fnames[i++] = "SubDomains";
	_assert_(i==numfields);

	/*Initialize Matlab structure*/
	dataref=mxCreateStructArray(ndim,dimensions,numfields,fnames);

	/*set each matlab each field*/
	i=0;
	i++; SetStructureField(dataref,"Vertices",        bamggeom->VerticesSize[0],        bamggeom->VerticesSize[1],        bamggeom->Vertices);
	i++; SetStructureField(dataref,"Edges",           bamggeom->EdgesSize[0],           bamggeom->EdgesSize[1],           bamggeom->Edges);
	i++; SetStructureField(dataref,"TangentAtEdges",  bamggeom->TangentAtEdgesSize[0],  bamggeom->TangentAtEdgesSize[1],  bamggeom->TangentAtEdges);
	i++; SetStructureField(dataref,"Corners",         bamggeom->CornersSize[0],         bamggeom->CornersSize[1],         bamggeom->Corners);
	i++; SetStructureField(dataref,"RequiredVertices",bamggeom->RequiredVerticesSize[0],bamggeom->RequiredVerticesSize[1],bamggeom->RequiredVertices);
	i++; SetStructureField(dataref,"RequiredEdges",   bamggeom->RequiredEdgesSize[0],   bamggeom->RequiredEdgesSize[1],   bamggeom->RequiredEdges);
	i++; SetStructureField(dataref,"CrackedEdges",    bamggeom->CrackedEdgesSize[0],    bamggeom->CrackedEdgesSize[1],    bamggeom->CrackedEdges);
	i++; SetStructureField(dataref,"SubDomains",      bamggeom->SubDomainsSize[0],      bamggeom->SubDomainsSize[1],      bamggeom->SubDomains);
	_assert_(i==numfields);

	/*Assign output*/
	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,BamgMesh* bamgmesh){{{*/
void WriteData(mxArray** pdataref,BamgMesh* bamgmesh){

	/*Intermediary*/
	int         i;
	mxArray    *dataref           = NULL;
	const int   numfields         = 16;
	const char *fnames[numfields];
	mwSize      ndim              = 2;
	mwSize      dimensions[2]     = {1,1};

	/*Initialize field names*/
	i=0;
	fnames[i++] = "Vertices";
	fnames[i++] = "Edges";
	fnames[i++] = "Triangles";
	fnames[i++] = "IssmEdges";
	fnames[i++] = "IssmSegments";
	fnames[i++] = "VerticesOnGeomVertex";
	fnames[i++] = "VerticesOnGeomEdge";
	fnames[i++] = "EdgesOnGeomEdge";
	fnames[i++] = "SubDomains";
	fnames[i++] = "SubDomainsFromGeom";
	fnames[i++] = "ElementConnectivity";
	fnames[i++] = "NodalConnectivity";
	fnames[i++] = "NodalElementConnectivity";
	fnames[i++] = "CrackedVertices";
	fnames[i++] = "CrackedEdges";
	fnames[i++] = "PreviousNumbering";
	_assert_(i==numfields);

	/*Initialize Matlab structure*/
	dataref=mxCreateStructArray(ndim,dimensions,numfields,fnames);

	/*set each matlab each field*/
	i=0;
	i++; SetStructureField(dataref,"Vertices",bamgmesh->VerticesSize[0], bamgmesh->VerticesSize[1],bamgmesh->Vertices);
	i++; SetStructureField(dataref,"Edges", bamgmesh->EdgesSize[0],bamgmesh->EdgesSize[1], bamgmesh->Edges);
	i++; SetStructureField(dataref,"Triangles", bamgmesh->TrianglesSize[0],bamgmesh->TrianglesSize[1], bamgmesh->Triangles);
	i++; SetStructureField(dataref,"IssmEdges", bamgmesh->IssmEdgesSize[0],bamgmesh->IssmEdgesSize[1], bamgmesh->IssmEdges);
	i++; SetStructureField(dataref,"IssmSegments",bamgmesh->IssmSegmentsSize[0], bamgmesh->IssmSegmentsSize[1],bamgmesh->IssmSegments);
	i++; SetStructureField(dataref,"VerticesOnGeomVertex",bamgmesh->VerticesOnGeomVertexSize[0],bamgmesh->VerticesOnGeomVertexSize[1], bamgmesh->VerticesOnGeomVertex);
	i++; SetStructureField(dataref,"VerticesOnGeomEdge",bamgmesh->VerticesOnGeomEdgeSize[0],bamgmesh->VerticesOnGeomEdgeSize[1], bamgmesh->VerticesOnGeomEdge);
	i++; SetStructureField(dataref,"EdgesOnGeomEdge", bamgmesh->EdgesOnGeomEdgeSize[0], bamgmesh->EdgesOnGeomEdgeSize[1],bamgmesh->EdgesOnGeomEdge);
	i++; SetStructureField(dataref,"SubDomains",bamgmesh->SubDomainsSize[0], bamgmesh->SubDomainsSize[1],bamgmesh->SubDomains);
	i++; SetStructureField(dataref,"SubDomainsFromGeom", bamgmesh->SubDomainsFromGeomSize[0], bamgmesh->SubDomainsFromGeomSize[1],bamgmesh->SubDomainsFromGeom);
	i++; SetStructureField(dataref,"ElementConnectivity",bamgmesh->ElementConnectivitySize[0],bamgmesh->ElementConnectivitySize[1], bamgmesh->ElementConnectivity);
	i++; SetStructureField(dataref,"NodalConnectivity",bamgmesh->NodalConnectivitySize[0],bamgmesh->NodalConnectivitySize[1], bamgmesh->NodalConnectivity);
	i++; SetStructureField(dataref,"NodalElementConnectivity", bamgmesh->NodalElementConnectivitySize[0], bamgmesh->NodalElementConnectivitySize[1],bamgmesh->NodalElementConnectivity);
	i++; SetStructureField(dataref,"CrackedVertices", bamgmesh->CrackedVerticesSize[0],bamgmesh->CrackedVerticesSize[1], bamgmesh->CrackedVertices);
	i++; SetStructureField(dataref,"CrackedEdges",bamgmesh->CrackedEdgesSize[0], bamgmesh->CrackedEdgesSize[1],bamgmesh->CrackedEdges);
	i++; SetStructureField(dataref,"PreviousNumbering",bamgmesh->VerticesSize[0],1,bamgmesh->PreviousNumbering);
	_assert_(i==numfields);

	/*Assign output*/
	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,IssmDenseMat<double>* matrix){{{*/
void WriteData(mxArray** pdataref,IssmDenseMat<double>* matrix){

	int      i,j;
	int      rows,cols;
	mxArray *dataref     = NULL;
	double  *matrix_ptr  = NULL;
	double  *tmatrix_ptr = NULL;

	if(matrix){

		matrix_ptr=matrix->ToSerial();
		matrix->GetSize(&rows,&cols);

		/*Now transpose the matrix and allocate with Matlab's memory manager: */
		tmatrix_ptr=(double*)mxMalloc(rows*cols*sizeof(double));
		for(i=0;i<rows;i++){
			for(j=0;j<cols;j++){
				tmatrix_ptr[j*rows+i]=matrix_ptr[i*cols+j];
			}
		}

		/*create matlab matrix: */
		dataref=mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,rows); 
		mxSetN(dataref,cols);
		mxSetPr(dataref,tmatrix_ptr);

		/*Free resources:*/
		xDelete<double>(matrix_ptr);
	}
	else{
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
	}

	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,IssmSeqVec<double>* vector){{{*/
void WriteData(mxArray** pdataref,IssmSeqVec<double>* vector){

	mxArray* dataref=NULL;
	double*  vector_ptr=NULL;
	double*  vector_matlab=NULL;
	int      rows;

	if(vector){
		/*call toolkit routine: */
		vector_ptr=vector->ToMPISerial();
		vector->GetSize(&rows);

		/*now create the matlab vector with Matlab's memory manager */
		vector_matlab=(double*)mxMalloc(rows*sizeof(double));
		for(int i=0;i<rows;i++) vector_matlab[i]=vector_ptr[i];

		dataref = mxCreateDoubleMatrix(0,0,mxREAL);                         
		mxSetM(dataref,rows);
		mxSetN(dataref,1);                                                                                          
		mxSetPr(dataref,vector_matlab);           
	}
	else{
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
	}

	/*Clean-up and return*/
	xDelete<double>(vector_ptr);
	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,RiftStruct* riftstruct){{{*/
void WriteData(mxArray** pdataref,RiftStruct* riftstruct){

	/*Intermediary*/
	int         i;
	mxArray    *dataref           = NULL;
	const int   numfields         = 10;
	const char *fnames[numfields];
	mwSize      ndim              = 2;
	mwSize      dimensions[2]     = {1,1};

	/*Initialize field names*/
	i=0;
	fnames[i++] = "numsegs";
	fnames[i++] = "segments";
	fnames[i++] = "pairs";
	fnames[i++] = "tips";
	fnames[i++] = "penaltypairs";
	fnames[i++] = "fill";
	fnames[i++] = "friction";
	fnames[i++] = "fraction";
	fnames[i++] = "fractionincrement";
	fnames[i++] = "state";
	_assert_(i==numfields);

	/*Initialize matlab structure of dimension numrifts*/
	dimensions[0]=riftstruct->numrifts;
	dataref=mxCreateStructArray(ndim,dimensions,numfields,fnames);

	/*set each matlab each field*/
	for(int i=0;i<riftstruct->numrifts;i++){
		SetStructureFieldi(dataref,i,"numsegs"          ,riftstruct->riftsnumsegments[i]);
		SetStructureFieldi(dataref,i,"segments"         ,riftstruct->riftsnumsegments[i]    ,3,riftstruct->riftssegments[i]);
		SetStructureFieldi(dataref,i,"pairs"            ,riftstruct->riftsnumpairs[i]       ,2,riftstruct->riftspairs[i]);
		SetStructureFieldi(dataref,i,"tips"             ,1                                  ,2,&riftstruct->riftstips[2*i]);
		SetStructureFieldi(dataref,i,"penaltypairs"     ,riftstruct->riftsnumpenaltypairs[i],7,riftstruct->riftspenaltypairs[i]);
		SetStructureFieldi(dataref,i,"fill"             ,"Ice");
		SetStructureFieldi(dataref,i,"friction"         ,0);
		SetStructureFieldi(dataref,i,"fraction"         ,0.);
		SetStructureFieldi(dataref,i,"fractionincrement",0.1);
		SetStructureFieldi(dataref,i,"state"            ,riftstruct->riftsnumpenaltypairs[i],1,riftstruct->state[i]);
	}

	/*Assign output*/
	*pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(mxArray** pdataref,Contours* contours){{{*/
void WriteData(mxArray** pdataref,Contours* contours){

	/*Intermediary*/

	int         i;
	mxArray    *dataref           = NULL;
	const int   numfields         = 6;
	const char *fnames[numfields];
	mwSize      ndim              = 2;
	mwSize      dimensions[2]     = {1,1};
	char        id[100];

	/*Initialize field names*/
	i=0;
	fnames[i++] = "name";
	fnames[i++] = "nods";
	fnames[i++] = "density";
	fnames[i++] = "x";
	fnames[i++] = "y";
	fnames[i++] = "closed";
	_assert_(i==numfields);

	/*Initialize matlab structure of dimension numrifts*/
	dimensions[0]=contours->Size();
	dataref=mxCreateStructArray(ndim,dimensions,numfields,fnames);

	/*set each matlab each field*/
	for(int i=0;i<contours->Size();i++){
		Contour<IssmPDouble>* contour=(Contour<IssmPDouble>*)contours->GetObjectByOffset(i);

		/*create a name for this contour from the contour id and set it: */
		sprintf(id,"%i",contour->id);
		SetStructureFieldi(dataref,i,"name",id);

		/*number of nods: */
		SetStructureFieldi(dataref,i,"nods"         ,contour->nods);

		/*density: */
		SetStructureFieldi(dataref,i,"density"            ,1);

		/*x and y: */
		SetStructureFieldi(dataref,i,"x"             ,contour->nods, 1,contour->x);
		SetStructureFieldi(dataref,i,"y"             ,contour->nods, 1,contour->y);

		/*closed: */
		SetStructureFieldi(dataref,i,"closed"         ,(int)contour->closed);
	}

	/*Assign output*/
	*pdataref=dataref;
}
/*}}}*/

/*Toolkit*/
/*FUNCTION SetStructureField(mxArray* dataref,const char* fieldname,int M,int N,double* fieldpointer){{{*/
void SetStructureField(mxArray* dataref,const char* fieldname,int M,int N,double* fieldpointer){

	mxArray* field = NULL;

	/*Convert field*/
	WriteData(&field,fieldpointer,M,N);

	/*Assign to structure*/
	mxSetField(dataref,0,fieldname,field);
}
/*}}}*/
/*FUNCTION SetStructureFieldi(mxArray* dataref,const char* fieldname,const char* string) {{{*/
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,const char* string){

	mxArray* field = NULL;

	/*Convert field*/
	WriteData(&field,string);

	/*Assign to structure*/
	mxSetField(dataref,i,fieldname,field);
}
/*}}}*/
/*FUNCTION SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int M,int N,double* fieldpointer){{{*/
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int M,int N,double* fieldpointer){

	mxArray* field = NULL;

	/*Convert field*/
	WriteData(&field,fieldpointer,M,N);

	/*Assign to structure*/
	mxSetField(dataref,i,fieldname,field);
}
/*}}}*/
/*FUNCTION SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int M,int N,int* fieldpointer){{{*/
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int M,int N,int* fieldpointer){

	mxArray* field = NULL;

	/*Convert field*/
	WriteData(&field,fieldpointer,M,N);

	/*Assign to structure*/
	mxSetField(dataref,i,fieldname,field);
}
/*}}}*/
/*FUNCTION SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int field){{{*/
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int fieldin){

	mxArray* field = NULL;

	/*Convert field*/
	WriteData(&field,fieldin);

	/*Assign to structure*/
	mxSetField(dataref,i,fieldname,field);
}
/*}}}*/
/*FUNCTION SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,double field){{{*/
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,double fieldin){

	mxArray* field = NULL;

	/*Convert field*/
	WriteData(&field,fieldin);

	/*Assign to structure*/
	mxSetField(dataref,i,fieldname,field);
}
/*}}}*/
