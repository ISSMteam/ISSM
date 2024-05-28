/* \file FetchJavascriptData.cpp:
 * \brief: general I/O interface to fetch data in javascript
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./javascriptio.h"
#include <cstring>

/*Primitive data types*/
/*FUNCTION FetchData(char** pstring, char* string){{{*/
void FetchData(char** pstring, char* stringin){
	char* string=NULL;

	string=xNew<char>(strlen(stringin)+1); xMemCpy<char>(string,stringin,strlen(stringin)+1);

	*pstring=string;
}
/*}}}*/
/*FUNCTION FetchData(double* pscalar, double scalar){{{*/
void FetchData(double* pscalar, double scalar){
	*pscalar=scalar;
}
/*}}}*/
/*FUNCTION FetchData(int* pinteger, int integer){{{*/
void FetchData(int* pinteger, int integer){
	*pinteger=integer;
}
/*}}}*/
/*FUNCTION FetchData(double **pvector, double* vectorin, int nods){{{*/
void FetchData(double** pvector, double* vectorin, int nods){
	double* vector=NULL;

	vector=xNew<IssmPDouble>(nods); xMemCpy<IssmPDouble>(vector,vectorin,nods);

	*pvector=vector;
}
/*}}}*/
/*FUNCTION FetchData(double **pvector, int* pnods, double* vectorin, int nods){{{*/
void FetchData(double** pvector, int* pnods, double* vectorin, int nods){
	double* vector=NULL;

	vector=xNew<IssmPDouble>(nods); xMemCpy<IssmPDouble>(vector,vectorin,nods);

	*pvector=vector;
	*pnods=nods;
}
/*}}}*/
/*FUNCTION FetchData(double **pmatrix, int* pM, int* matrix, int M, int N){{{*/
void FetchData(double **pmatrix, int* pM, int* matrixin, int M, int N){
	double*  outmatrix=NULL;
	int      outmatrix_rows;

	if(M == 0 || N == 0){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix=NULL;
	}
	else if (pmatrix && matrixin){
		outmatrix_rows=M;
		outmatrix=xNew<IssmPDouble>(M*N);
		for(int i=0;i<M*N;i++){outmatrix[i]=(IssmPDouble)matrixin[i];}
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM){*pM=outmatrix_rows;}
}
/*}}}*/
/*FUNCTION FetchData(double **pmatrix, int* pM, int* pN, int* matrix, int M, int N){{{*/
void FetchData(double **pmatrix, int* pM, int* pN, int* matrixin, int M, int N){
	double*  outmatrix=NULL;
	int      outmatrix_rows, outmatrix_cols;

	if(M == 0 || N == 0){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if (pmatrix && matrixin){
		outmatrix_rows=M;
		outmatrix_cols=N;
		outmatrix=xNew<IssmPDouble>(M*N);
		for(int i=0;i<M*N;i++){outmatrix[i]=(IssmPDouble)matrixin[i];}
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM){*pM=outmatrix_rows;}
	if (pN){*pN=outmatrix_cols;}
}
/*}}}*/
/*FUNCTION FetchData(double **pmatrix, int* pM, double* matrix, int M, int N){{{*/
void FetchData(double **pmatrix, int* pM, double* matrixin, int M, int N){
	double*  outmatrix=NULL;
	int      outmatrix_rows;

	if(M == 0 || N == 0){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix=NULL;
	}
	else if (pmatrix && matrixin){
		outmatrix_rows=M;
		outmatrix=xNew<IssmPDouble>(M*N); xMemCpy<IssmPDouble>(outmatrix,matrixin,M*N);
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM){*pM=outmatrix_rows;}
}
/*}}}*/
/*FUNCTION FetchData(double **pmatrix, int* pM, int* pN, double* matrix, int M, int N){{{*/
void FetchData(double **pmatrix, int* pM, int* pN, double* matrixin, int M, int N){
	double*  outmatrix=NULL;
	int      outmatrix_rows, outmatrix_cols;

	if(M == 0 || N == 0){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if (pmatrix && matrixin){
		outmatrix_rows=M;
		outmatrix_cols=N;
		outmatrix=xNew<IssmPDouble>(M*N); xMemCpy<IssmPDouble>(outmatrix,matrixin,M*N);
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM){*pM=outmatrix_rows;}
	if (pN){*pN=outmatrix_cols;}
}
/*}}}*/
/*FUNCTION FetchData(int **pmatrix, int* pM, int* matrix, int M, int N){{{*/
void FetchData(int **pmatrix, int* pM, int* matrixin, int M, int N){
	int*     outmatrix=NULL;
	int      outmatrix_rows;

	if(M == 0 || N == 0){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix=NULL;
	}
	else if (pmatrix && matrixin){
		outmatrix_rows=M;
		outmatrix=xNew<int>(M*N); xMemCpy<int>(outmatrix,matrixin,M*N);
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM){*pM=outmatrix_rows;}
}
/*}}}*/
/*FUNCTION FetchData(int **pmatrix, int* pM, int* pN, int* matrix, int M, int N){{{*/
void FetchData(int **pmatrix, int* pM, int* pN, int* matrixin, int M, int N){
	int*     outmatrix=NULL;
	int      outmatrix_rows, outmatrix_cols;

	if(M == 0 || N == 0){
		/*Nothing to pick up. Just initialize matrix pointer to NULL: */
		outmatrix_rows=0;
		outmatrix_cols=0;
		outmatrix=NULL;
	}
	else if (pmatrix && matrixin){
		outmatrix_rows=M;
		outmatrix_cols=N;
		outmatrix=xNew<int>(M*N); xMemCpy<int>(outmatrix,matrixin,M*N);
	}

	/*Assign output pointers:*/
	*pmatrix=outmatrix;
	if (pM){*pM=outmatrix_rows;}
	if (pN){*pN=outmatrix_cols;}
}
/*}}}*/
/*ISSM objects*/
/*FUNCTION FetchData(Contours** pcontours, double* x, double* y, int nods){{{*/
void FetchData(Contours** pcontours, double* x, double* y, int nods){
	int             numcontours, index, test1, test2;
	char            *contourname = NULL;
	Contours        *contours    = NULL;
	Contour<double> *contouri    = NULL;

	/*only 1 contour for now: */
	contours=new Contours();

	if (nods){
		contouri=new Contour<double>();
		contouri->nods=nods;
		contouri->x=xNew<IssmPDouble>(nods); xMemCpy<IssmPDouble>(contouri->x,x,nods);
		contouri->y=xNew<IssmPDouble>(nods); xMemCpy<IssmPDouble>(contouri->y,y,nods);
		contours->AddObject(contouri);
	}

	*pcontours=contours;
}
/*}}}*/
/*FUNCTION FetchData(BamgGeom** pbamggeom, int* VerticesSize, double* Vertices, int* EdgesSize, double* Edges, int* CornersSize, double* Corners, int* RequiredVerticesSize, double* RequiredVertices, int* RequiredEdgesSize, double* RequiredEdges, int* CrackedEdgesSize, double* CrackedEdges, int* SubDomainsSize, double* SubDomains){{{*/
void FetchData(BamgGeom** pbamggeom, int* VerticesSize, double* Vertices, int* EdgesSize, double* Edges, int* CornersSize, double* Corners, int* RequiredVerticesSize, double* RequiredVertices, int* RequiredEdgesSize, double* RequiredEdges, int* CrackedEdgesSize, double* CrackedEdges, int* SubDomainsSize, double* SubDomains){

	/*Initialize output*/
	BamgGeom* bamggeom=new BamgGeom();

	/*Fetch all fields*/
	FetchData(&bamggeom->Vertices, &bamggeom->VerticesSize[0], &bamggeom->VerticesSize[1], Vertices, VerticesSize[0], VerticesSize[1]);
	FetchData(&bamggeom->Edges, &bamggeom->EdgesSize[0], &bamggeom->EdgesSize[1], Edges, EdgesSize[0], EdgesSize[1]);
	FetchData(&bamggeom->Corners, &bamggeom->CornersSize[0], &bamggeom->CornersSize[1], Corners, CornersSize[0], CornersSize[1]);
	FetchData(&bamggeom->RequiredVertices, &bamggeom->RequiredVerticesSize[0], &bamggeom->RequiredVerticesSize[1], RequiredVertices, RequiredVerticesSize[0], RequiredVerticesSize[1]);
	FetchData(&bamggeom->RequiredEdges, &bamggeom->RequiredEdgesSize[0], &bamggeom->RequiredEdgesSize[1], RequiredEdges, RequiredEdgesSize[0], RequiredEdgesSize[1]);
	FetchData(&bamggeom->CrackedEdges, &bamggeom->CrackedEdgesSize[0], &bamggeom->CrackedEdgesSize[1], CrackedEdges, CrackedEdgesSize[0], CrackedEdgesSize[1]);
	FetchData(&bamggeom->SubDomains, &bamggeom->SubDomainsSize[0], &bamggeom->SubDomainsSize[1], SubDomains, SubDomainsSize[0], SubDomainsSize[1]);

	/*Assign output pointers:*/
	*pbamggeom=bamggeom;
}
/*}}}*/
/*FUNCTION FetchData(BamgMesh** pbamgmesh, int* VerticesSize, double* Vertices, int* EdgesSize, double* Edges, int* TrianglesSize, double* Triangles, int* CrackedEdgesSize, double* CrackedEdges, int* VerticesOnGeomEdgeSize, double* VerticesOnGeomEdge, int* VerticesOnGeomVertexSize, double* VerticesOnGeomVertex, int* EdgesOnGeomEdgeSize, double* EdgesOnGeomEdge, int* IssmSegmentsSize, double* IssmSegments){{{*/
void FetchData(BamgMesh** pbamgmesh, int* VerticesSize, double* Vertices, int* EdgesSize, double* Edges, int* TrianglesSize, double* Triangles, int* CrackedEdgesSize, double* CrackedEdges, int* VerticesOnGeomEdgeSize, double* VerticesOnGeomEdge, int* VerticesOnGeomVertexSize, double* VerticesOnGeomVertex, int* EdgesOnGeomEdgeSize, double* EdgesOnGeomEdge, int* IssmSegmentsSize, double* IssmSegments){

	/*Initialize output*/
	BamgMesh* bamgmesh=new BamgMesh();

	/*Fetch all fields*/
	FetchData(&bamgmesh->Vertices, &bamgmesh->VerticesSize[0], &bamgmesh->VerticesSize[1], Vertices, VerticesSize[0], VerticesSize[1]);
	FetchData(&bamgmesh->Edges, &bamgmesh->EdgesSize[0], &bamgmesh->EdgesSize[1], Edges, EdgesSize[0], EdgesSize[1]);
	FetchData(&bamgmesh->Triangles, &bamgmesh->TrianglesSize[0], &bamgmesh->TrianglesSize[1], Triangles, TrianglesSize[0], TrianglesSize[1]);
	FetchData(&bamgmesh->CrackedEdges, &bamgmesh->CrackedEdgesSize[0], &bamgmesh->CrackedEdgesSize[1], CrackedEdges, CrackedEdgesSize[0], CrackedEdgesSize[1]);
	FetchData(&bamgmesh->VerticesOnGeomEdge, &bamgmesh->VerticesOnGeomEdgeSize[0], &bamgmesh->VerticesOnGeomEdgeSize[1], VerticesOnGeomEdge, VerticesOnGeomEdgeSize[0], VerticesOnGeomEdgeSize[1]);
	FetchData(&bamgmesh->VerticesOnGeomVertex, &bamgmesh->VerticesOnGeomVertexSize[0], &bamgmesh->VerticesOnGeomVertexSize[1], VerticesOnGeomVertex, VerticesOnGeomVertexSize[0], VerticesOnGeomVertexSize[1]);
	FetchData(&bamgmesh->EdgesOnGeomEdge, &bamgmesh->EdgesOnGeomEdgeSize[0], &bamgmesh->EdgesOnGeomEdgeSize[1], EdgesOnGeomEdge, EdgesOnGeomEdgeSize[0], EdgesOnGeomEdgeSize[1]);
	FetchData(&bamgmesh->IssmSegments, &bamgmesh->IssmSegmentsSize[0], &bamgmesh->IssmSegmentsSize[1], IssmSegments, IssmSegmentsSize[0], IssmSegmentsSize[1]);

	/*Assign output pointers:*/
	*pbamgmesh=bamgmesh;
}
/*}}}*/
/*FUNCTION FetchData(BamgOpts** pbamgopts, double anisomax, double cutoff, double coeff, double errg, double gradation, int Hessiantype, int maxnbv, double maxsubdiv, int Metrictype, int nbjacobi, int nbsmooth, double omega, double power, int verbose, int Crack, int KeepVertices, int splitcorners, double hmin, double hmax, int* hminVerticesSize, double* hminVertices, int* hmaxVerticesSize, double* hmaxVertices, int hVerticesLength, double* hVertices, int* metricSize, double* metric, int* fieldSize, double* field, int* errSize, double* err){{{*/
void FetchData(BamgOpts** pbamgopts, double anisomax, double coeff, double cutoff, double errg, double gradation, int Hessiantype, int maxnbv, double maxsubdiv, int Metrictype, int nbjacobi, int nbsmooth, double omega, double power, int verbose, int Crack, int KeepVertices, int splitcorners, double hmin, double hmax, int* hminVerticesSize, double* hminVertices, int* hmaxVerticesSize, double* hmaxVertices, int hVerticesLength, double* hVertices, int* metricSize, double* metric, int* fieldSize, double* field, int* errSize, double* err){

	BamgOpts *bamgopts      = new BamgOpts();

	/*Parameters*/
	bamgopts->anisomax	    = anisomax;
	bamgopts->coeff	        = coeff;
	bamgopts->cutoff    	= cutoff;
	bamgopts->errg	        = errg;
	bamgopts->gradation	    = gradation;
	bamgopts->Hessiantype	= Hessiantype;
	bamgopts->maxnbv	    = maxnbv;
	bamgopts->maxsubdiv	    = maxsubdiv;
	bamgopts->Metrictype	= Metrictype;
	bamgopts->nbjacobi	    = nbjacobi;
	bamgopts->nbsmooth	    = nbsmooth;
	bamgopts->omega	        = omega;
	bamgopts->power	        = power;
	bamgopts->verbose	    = verbose;

	/*Flags*/
	bamgopts->Crack	        = Crack;
	bamgopts->KeepVertices	= KeepVertices;
	bamgopts->splitcorners	= splitcorners;

	/*Metric related*/
	bamgopts->hmin	        = hmin;
	bamgopts->hmax       	= hmax;
	FetchData(&bamgopts->hminVertices, &bamgopts->hminVerticesSize[0], &bamgopts->hminVerticesSize[1], hminVertices, hminVerticesSize[0], hminVerticesSize[1]);
	FetchData(&bamgopts->hmaxVertices, &bamgopts->hmaxVerticesSize[0], &bamgopts->hmaxVerticesSize[1], hmaxVertices, hmaxVerticesSize[0], hmaxVerticesSize[1]);
	FetchData(&bamgopts->hVertices, &bamgopts->hVerticesLength, hVertices, hVerticesLength);
	FetchData(&bamgopts->field, &bamgopts->fieldSize[0], &bamgopts->fieldSize[1], field, fieldSize[0], fieldSize[1]);
	FetchData(&bamgopts->metric, &bamgopts->metricSize[0], &bamgopts->metricSize[1], metric, metricSize[0], metricSize[1]);
	FetchData(&bamgopts->err, &bamgopts->errSize[0], &bamgopts->errSize[1], err, errSize[0], errSize[1]);

	/*Additional checks*/
	bamgopts->Check();

	/*Assign output pointers:*/
	*pbamgopts=bamgopts;
}
/*}}}*/
/*FUNCTION FetchData(Options** poptions, int NRHS, int nrhs, const char* optionname, double optionvalue){{{*/
void FetchData(Options** poptions, int NRHS, int nrhs, const char* optionname, double optionvalue){

	/*Initialize output*/
	Options* options=new Options();

	/*check and parse the name  */
	GenericOption<double> *odouble=new GenericOption<double>();
	odouble=new GenericOption<double>();
	odouble->name=xNew<char>(strlen(optionname)+1);
	memcpy(odouble->name,optionname,(strlen(optionname)+1)*sizeof(char));
	odouble->value=optionvalue;
	odouble->size[0]=1;
	odouble->size[1]=1;
	options->AddOption((Option*)odouble);

	/*Assign output pointers:*/
	*poptions=options;
}
/*}}}*/
