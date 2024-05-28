/*\file javascriptio.h
 *s\brief: I/O for ISSM in javascript mode
 */

#ifndef _JAVASCRIPT_IO_H_
#define _JAVASCRIPT_IO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif 

#include "../include/javascriptincludes.h"
#include "../../../c/bamg/bamgobjects.h"
#include "../../../c/classes/classes.h"
#include "../../../c/toolkits/toolkits.h"
#include "../../../c/shared/shared.h"

void WriteData(IssmPDouble** pmatrix, int* pnel, int* matrix, int M,int N);
void WriteData(IssmPDouble** pmatrix, int* pM, int* pN, int* matrix, int M, int N);
void WriteData(IssmPDouble** pmatrix, int* pM, int* pN, IssmPDouble* matrix, int M, int N);
void WriteData(IssmPDouble** pmatrix, int** pSize, IssmPDouble* matrix, int* size);
void WriteData(IssmPDouble** px, int* pnods, int* vector, int M);
void WriteData(IssmPDouble** px, int* pnods, double* vector, int M);
void WriteData(char** pstring, char* stringin);
void WriteData(IssmPDouble* pdouble, IssmPDouble doublein);
void WriteData(IssmPDouble** pdataref, IssmSeqVec<double>* vector);
void WriteData(IssmPDouble** pdouble, void*);
void WriteData(int** VerticesSize, double** Vertices, int** EdgesSize, double** Edges, int** CornersSize, double** Corners, int** RequiredVerticesSize, double** RequiredVertices, int** RequiredEdgesSize, double** RequiredEdges, int** CrackedEdgesSize, double** CrackedEdges, int** SubDomainsSize, double** SubDomains, BamgGeom* bamggeom);
void WriteData(int** VerticesSize, double** Vertices, int** EdgesSize, double** Edges, int** TrianglesSize, double** Triangles, int** IssmEdgesSize, double** IssmEdges, int** IssmSegmentsSize, double** IssmSegments, int** VerticesOnGeomVertexSize, double** VerticesOnGeomVertex, int** VerticesOnGeomEdgeSize, double** VerticesOnGeomEdge, int** EdgesOnGeomEdgeSize, double** EdgesOnGeomEdge, int** SubDomainsSize, double** SubDomains, int** SubDomainsFromGeomSize, double** SubDomainsFromGeom, int** ElementConnectivitySize, double** ElementConnectivity, int** NodalConnectivitySize, double** NodalConnectivity, int** NodalElementConnectivitySize, double** NodalElementConnectivity, int** CrackedVerticesSize, double** CrackedVertices, int** CrackedEdgesSize, double** CrackedEdges, int** PreviousNumberingSize, double** PreviousNumbering, BamgMesh* bamgmesh);

void FetchData(char** pstring, char* stringin);
void FetchData(double* pscalar, double scalar);
void FetchData(int* pinteger,int integer);
void FetchData(double** pvector, double* vectorin, int nods);
void FetchData(double** pvector, int* pnods, double* vectorin, int nods);
void FetchData(double** pmatrix, int* pM, int* matrixin, int M, int N);
void FetchData(double** pmatrix, int* pM, int* pN, int* matrixin, int M, int N);
void FetchData(double** pmatrix, int* pM, double* matrixin, int M, int N);
void FetchData(double** pmatrix, int* pM, int* pN, double* matrixin, int M, int N);
void FetchData(int** pmatrix, int* pM, int* matrixin, int M, int N);
void FetchData(int** pmatrix, int* pM, int* pN, int* matrixin, int M, int N);
void FetchData(Contours** pcontours,double* x, double* y, int nods);
void FetchData(BamgGeom** pbamggeom, int* VerticesSize, double* Vertices, int* EdgesSize, double* Edges, int* CornersSize, double* Corners, int* RequiredVerticesSize, double* RequiredVertices, int* RequiredEdgesSize, double* RequiredEdges, int* CrackedEdgesSize, double* CrackedEdges, int* SubDomainsSize, double* SubDomains);
void FetchData(BamgMesh** pbamgmesh, int* VerticesSize, double* Vertices, int* EdgesSize, double* Edges, int* TrianglesSize, double* Triangles, int* CrackedEdgesSize, double* CrackedEdges, int* VerticesOnGeomEdgeSize, double* VerticesOnGeomEdge, int* VerticesOnGeomVertexSize, double* VerticesOnGeomVertex, int* EdgesOnGeomEdgeSize, double* EdgesOnGeomEdge, int* IssmSegmentsSize, double* IssmSegments);
void FetchData(BamgOpts** pbamgopts, double anisomax, double cutoff, double coeff, double errg, double gradation, int Hessiantype, int maxnbv, double maxsubdiv, int Metrictype, int nbjacobi, int nbsmooth, double omega, double power, int verbose, int Crack, int KeepVertices, int splitcorners, double hmin, double hmax, int* hminVerticesSize, double* hminVertices, int* hmaxVerticesSize, double* hmaxVertices, int hVerticesLength, double* hVertices, int* metricSize, double* metric, int* fieldSize, double* field, int* errSize, double* err);
void FetchData(Options** poptions,int NRHS, int nrhs, const char* optionname, double optionvalue);
void FetchData(int* pinteger,int integer);

/*Print*/
void ApiPrintf(const char* string);
#endif	/* _IO_H_ */
