/* \file WriteJavascriptData.cpp:
 * \brief: general I/O interface to fetch data in javascript
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./javascriptio.h"
#include "./../../../c/datastructures/datastructures.h"

/*Primitive data types*/
/*FUNCTION WriteData(IssmPDouble** pmatrix, int* pnel, int* matrix, int M,int N){{{*/
void WriteData(IssmPDouble** pmatrix, int* pnel, int* matrix, int M,int N){

	if(pmatrix && matrix){

		/*Copy matrix: */
		IssmPDouble* dmatrix = xNew<IssmPDouble>(M*N); 
		for (int i=0;i<M*N;i++)dmatrix[i]=(IssmPDouble)matrix[i];
		*pmatrix=dmatrix;
		*pnel=M;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** pmatrix, int* pM, int* pN, , int* matrix, int M,int N){{{*/
void WriteData(IssmPDouble** pmatrix, int* pM, int* pN, int* matrix, int M, int N){

	if(pmatrix && matrix){

		/*Copy matrix: */
	    IssmPDouble* dmatrix = xNew<IssmPDouble>(M*N); 
		for (int i=0;i<M*N;i++) dmatrix[i]=(IssmPDouble)matrix[i];
		*pmatrix=dmatrix;
		*pM=M;
		*pN=N;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** pmatrix, int* pM, IssmPDouble* pN, , int* matrix, int M,int N){{{*/
void WriteData(IssmPDouble** pmatrix, int* pM, int* pN, IssmPDouble* matrix, int M, int N){

	if(pmatrix && matrix){

		/*Copy matrix: */
		IssmPDouble* dmatrix = xNew<IssmPDouble>(M*N); 
		for (int i=0;i<M*N;i++) dmatrix[i]=matrix[i];
		*pmatrix=dmatrix;
		*pM=M;
		*pN=N;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** pmatrix, int** pSize, IssmPDouble* matrix, int* size){{{*/
void WriteData(IssmPDouble** pmatrix, int** pSize, IssmPDouble* matrix, int* size){

    int M = size[0];
    int N = size[1];
    int* imatrix = xNew<int>(2);
    IssmPDouble* dmatrix = xNew<IssmPDouble>(M*N); 

    /*Copy matrix: */
    for (int i=0;i<2;i++) imatrix[i]=size[i];
    for (int i=0;i<M*N;i++) dmatrix[i]=matrix[i];
    *pmatrix=dmatrix;
    *pSize=imatrix;
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** px, int* pnods, double* vector, int M){{{*/
void WriteData(IssmPDouble** px, int* pnods, double* vector, int M){

	if(px && vector){

		IssmPDouble* dx=xNew<IssmPDouble>(M); 
		for(int i=0;i<M;i++)dx[i]=vector[i];
		*px=dx;
		*pnods=M;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** px, int* pnods, int* vector, int M){{{*/
void WriteData(IssmPDouble** px, int* pnods, int* vector, int M){

	if(px && vector){

		IssmPDouble* dx=xNew<IssmPDouble>(M); 
		for(int i=0;i<M;i++)dx[i]=(IssmPDouble)vector[i];
		*px=dx;
		*pnods=M;
	}
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble* pdouble, IssmSeqVec<double> vector){{{*/
void WriteData(IssmPDouble** pdataref, IssmSeqVec<double>* vector){

    double*  dataref=NULL;
    double*  vector_ptr=NULL;
    int      rows;

    if(vector){
        /*call toolkit routine: */
        vector_ptr=vector->ToMPISerial();
        vector->GetSize(&rows);

        /*now create the js vector */
		dataref=xNew<double>(rows); 
        for(int i=0;i<rows;i++) dataref[i]=vector_ptr[i];
    }

    /*Clean-up and return*/
    xDelete<double>(vector_ptr);
    *pdataref=dataref;
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble* pdouble, IssmPDouble double){{{*/
void WriteData(IssmPDouble* pdouble, IssmPDouble doublein){

	*pdouble=doublein;
}
/*}}}*/
/*FUNCTION WriteData(IssmPDouble** pdouble, void* nullptr){{{*/
void WriteData(IssmPDouble** pdouble, void*){
	//do nothing
}
/*}}}*/
/*FUNCTION WriteData(char** pstring, char* string){{{*/
void WriteData(char** pstring, char* stringin){

	char* string=xNew<char>(strlen(stringin)+1);
	xMemCpy<char>(string,stringin,strlen(stringin)+1);

	*pstring=string;
}
/*}}}*/

/*ISSM objects*/
/*FUNCTION WriteData(int** VerticesSize, double** Vertices, int** EdgesSize, double** Edges, int** CornersSize, double** Corners, int** RequiredVerticesSize, double** RequiredVertices, int** RequiredEdgesSize, double** RequiredEdges, int** CrackedEdgesSize, double** CrackedEdges, int** SubDomainsSize, double** SubDomains, BamgGeom* bamggeom){{{*/
void WriteData(int** VerticesSize, double** Vertices, int** EdgesSize, double** Edges, int** CornersSize, double** Corners, int** RequiredVerticesSize, double** RequiredVertices, int** RequiredEdgesSize, double** RequiredEdges, int** CrackedEdgesSize, double** CrackedEdges, int** SubDomainsSize, double** SubDomains, BamgGeom* bamggeom){

	/*Assign each field to output*/
    WriteData(Vertices, VerticesSize, bamggeom->Vertices, bamggeom->VerticesSize);
    WriteData(Edges, EdgesSize, bamggeom->Edges, bamggeom->EdgesSize);
    WriteData(Corners, CornersSize, bamggeom->Corners, bamggeom->CornersSize);
    WriteData(RequiredVertices, RequiredVerticesSize, bamggeom->RequiredVertices, bamggeom->RequiredVerticesSize);
    WriteData(RequiredEdges, RequiredEdgesSize, bamggeom->RequiredEdges, bamggeom->RequiredEdgesSize);
    WriteData(CrackedEdges, CrackedEdgesSize, bamggeom->CrackedEdges, bamggeom->CrackedEdgesSize);
    WriteData(SubDomains, SubDomainsSize, bamggeom->SubDomains, bamggeom->SubDomainsSize);
}
/*}}}*/
/*FUNCTION WriteData(int** VerticesSize, double** Vertices, int** EdgesSize, double** Edges, int** TrianglesSize, double** Triangles, int** IssmEdgesSize, double** IssmEdges, int** IssmSegmentsSize, double** IssmSegments, int** VerticesOnGeomVertexSize, double** VerticesOnGeomVertex, int** VerticesOnGeomEdgeSize, double** VerticesOnGeomEdge, int** EdgesOnGeomEdgeSize, double** EdgesOnGeomEdge, int** SubDomainsSize, double** SubDomains, int** SubDomainsFromGeomSize, double** SubDomainsFromGeom, int** ElementConnectivitySize, double** ElementConnectivity, int** NodalConnectivitySize, double** NodalConnectivity, int** NodalElementConnectivitySize, double** NodalElementConnectivity, int** CrackedVerticesSize, double** CrackedVertices, int** CrackedEdgesSize, double** CrackedEdges, int** PreviousNumberingSize, double** PreviousNumbering, BamgMesh* bamgmesh){{{*/
void WriteData(int** VerticesSize, double** Vertices, int** EdgesSize, double** Edges, int** TrianglesSize, double** Triangles, int** IssmEdgesSize, double** IssmEdges, int** IssmSegmentsSize, double** IssmSegments, int** VerticesOnGeomVertexSize, double** VerticesOnGeomVertex, int** VerticesOnGeomEdgeSize, double** VerticesOnGeomEdge, int** EdgesOnGeomEdgeSize, double** EdgesOnGeomEdge, int** SubDomainsSize, double** SubDomains, int** SubDomainsFromGeomSize, double** SubDomainsFromGeom, int** ElementConnectivitySize, double** ElementConnectivity, int** NodalConnectivitySize, double** NodalConnectivity, int** NodalElementConnectivitySize, double** NodalElementConnectivity, int** CrackedVerticesSize, double** CrackedVertices, int** CrackedEdgesSize, double** CrackedEdges, int** PreviousNumberingSize, double** PreviousNumbering, BamgMesh* bamgmesh){

	/*Assign each field to output*/
    WriteData(Vertices, VerticesSize, bamgmesh->Vertices, bamgmesh->VerticesSize);
    WriteData(Edges, EdgesSize, bamgmesh->Edges, bamgmesh->EdgesSize);
    WriteData(Triangles, TrianglesSize, bamgmesh->Triangles, bamgmesh->TrianglesSize);
    WriteData(IssmEdges, IssmEdgesSize, bamgmesh->IssmEdges, bamgmesh->IssmEdgesSize);
    WriteData(IssmSegments, IssmSegmentsSize, bamgmesh->IssmSegments, bamgmesh->IssmSegmentsSize);
    WriteData(VerticesOnGeomVertex, VerticesOnGeomVertexSize, bamgmesh->VerticesOnGeomVertex, bamgmesh->VerticesOnGeomVertexSize);
    WriteData(VerticesOnGeomEdge, VerticesOnGeomEdgeSize, bamgmesh->VerticesOnGeomEdge, bamgmesh->VerticesOnGeomEdgeSize);
    WriteData(EdgesOnGeomEdge, EdgesOnGeomEdgeSize, bamgmesh->EdgesOnGeomEdge, bamgmesh->EdgesOnGeomEdgeSize);
    WriteData(SubDomains, SubDomainsSize, bamgmesh->SubDomains, bamgmesh->SubDomainsSize);
    WriteData(SubDomainsFromGeom, SubDomainsFromGeomSize, bamgmesh->SubDomainsFromGeom, bamgmesh->SubDomainsFromGeomSize);
    WriteData(ElementConnectivity, ElementConnectivitySize, bamgmesh->ElementConnectivity, bamgmesh->ElementConnectivitySize);
    WriteData(NodalConnectivity, NodalConnectivitySize, bamgmesh->NodalConnectivity, bamgmesh->NodalConnectivitySize);
    WriteData(NodalElementConnectivity, NodalElementConnectivitySize, bamgmesh->NodalElementConnectivity, bamgmesh->NodalElementConnectivitySize);
    WriteData(CrackedVertices, CrackedVerticesSize, bamgmesh->CrackedVertices, bamgmesh->CrackedVerticesSize);
    WriteData(CrackedEdges, CrackedEdgesSize, bamgmesh->CrackedEdges, bamgmesh->CrackedEdgesSize);
    WriteData(PreviousNumbering, PreviousNumberingSize, bamgmesh->PreviousNumbering, bamgmesh->VerticesSize); //PreviousNumbering just resuses Vertices' Size
}
/*}}}*/
