/*
	BamgMesher.h
*/

#ifndef _BAMG_MESHER_H_
#define _BAMG_MESHER_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
	#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*For python modules: needs to come before header files inclusion*/
#ifdef _HAVE_PYTHON_
#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#endif

#include "../../c/main/globals.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"
#include "../bindings.h"

#undef __FUNCT__
#define __FUNCT__  "BamgMesher"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define BAMGMESHIN  prhs[0]
#define BAMGGEOMIN  prhs[1]
#define BAMGOPTIONS prhs[2]
/* serial output macros: */
#define BAMGMESHOUT (mxArray**)&plhs[0]
#define BAMGGEOMOUT (mxArray**)&plhs[1]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define BAMGMESHIN  PyTuple_GetItem(args,0)
#define BAMGGEOMIN  PyTuple_GetItem(args,1)
#define BAMGOPTIONS PyTuple_GetItem(args,2)
/* serial output macros: */
#define BAMGMESHOUT output,0
#define BAMGGEOMOUT output,1
#endif

#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define BAMGMESHIN VerticesSize_mesh_in, Vertices_mesh_in, EdgesSize_mesh_in, Edges_mesh_in, TrianglesSize_mesh_in, Triangles_mesh_in, CrackedEdgesSize_mesh_in, CrackedEdges_mesh_in, VerticesOnGeomEdgeSize_mesh_in, VerticesOnGeomEdge_mesh_in, VerticesOnGeomVertexSize_mesh_in, VerticesOnGeomVertex_mesh_in, EdgesOnGeomEdgeSize_mesh_in, EdgesOnGeomEdge_mesh_in, IssmSegmentsSize_mesh_in, IssmSegments_mesh_in
#define BAMGGEOMIN VerticesSize_geom_in, Vertices_geom_in, EdgesSize_geom_in, Edges_geom_in, CornersSize_geom_in, Corners_geom_in, RequiredVerticesSize_geom_in, RequiredVertices_geom_in, RequiredEdgesSize_geom_in, RequiredEdges_geom_in, CrackedEdgesSize_geom_in, CrackedEdges_geom_in, SubDomainsSize_geom_in, SubDomains_geom_in
#define BAMGOPTIONS anisomax, cutoff, coeff, errg, gradation, Hessiantype, maxnbv, maxsubdiv, Metrictype, nbjacobi, nbsmooth, omega, power, verbose, Crack, KeepVertices, splitcorners, hmin, hmax, hminVerticesSize, hminVertices, hmaxVerticesSize, hmaxVertices, hVerticesLength, hVertices, metricSize, metric, fieldSize, field, errSize, err
/* serial output macros: */
#define BAMGMESHOUT VerticesSize_mesh_out, Vertices_mesh_out, EdgesSize_mesh_out, Edges_mesh_out, TrianglesSize_mesh_out, Triangles_mesh_out, IssmEdgesSize_mesh_out, IssmEdges_mesh_out, IssmSegmentsSize_mesh_out, IssmSegments_mesh_out, VerticesOnGeomVertexSize_mesh_out, VerticesOnGeomVertex_mesh_out, VerticesOnGeomEdgeSize_mesh_out, VerticesOnGeomEdge_mesh_out, EdgesOnGeomEdgeSize_mesh_out, EdgesOnGeomEdge_mesh_out, SubDomainsSize_mesh_out, SubDomains_mesh_out, SubDomainsFromGeomSize_mesh_out, SubDomainsFromGeom_mesh_out, ElementConnectivitySize_mesh_out, ElementConnectivity_mesh_out, NodalConnectivitySize_mesh_out, NodalConnectivity_mesh_out, NodalElementConnectivitySize_mesh_out, NodalElementConnectivity_mesh_out, CrackedVerticesSize_mesh_out, CrackedVertices_mesh_out, CrackedEdgesSize_mesh_out, CrackedEdges_mesh_out, PreviousNumberingSize_mesh_out, PreviousNumbering_mesh_out
#define BAMGGEOMOUT VerticesSize_geom_out, Vertices_geom_out, EdgesSize_geom_out, Edges_geom_out, CornersSize_geom_out, Corners_geom_out, RequiredVerticesSize_geom_out, RequiredVertices_geom_out, RequiredEdgesSize_geom_out, RequiredEdges_geom_out, CrackedEdgesSize_geom_out, CrackedEdges_geom_out, SubDomainsSize_geom_out, SubDomains_geom_out
#define WRAPPER(modulename) extern "C" { int  BamgMesherModule(int** VerticesSize_mesh_out, double** Vertices_mesh_out, int** EdgesSize_mesh_out, double** Edges_mesh_out, int** TrianglesSize_mesh_out, double** Triangles_mesh_out, int** IssmEdgesSize_mesh_out, double** IssmEdges_mesh_out, int** IssmSegmentsSize_mesh_out, double** IssmSegments_mesh_out, int** VerticesOnGeomVertexSize_mesh_out, double** VerticesOnGeomVertex_mesh_out, int** VerticesOnGeomEdgeSize_mesh_out, double** VerticesOnGeomEdge_mesh_out, int** EdgesOnGeomEdgeSize_mesh_out, double** EdgesOnGeomEdge_mesh_out, int** SubDomainsSize_mesh_out, double** SubDomains_mesh_out, int** SubDomainsFromGeomSize_mesh_out, double** SubDomainsFromGeom_mesh_out, int** ElementConnectivitySize_mesh_out, double** ElementConnectivity_mesh_out, int** NodalConnectivitySize_mesh_out, double** NodalConnectivity_mesh_out, int** NodalElementConnectivitySize_mesh_out, double** NodalElementConnectivity_mesh_out, int** CrackedVerticesSize_mesh_out, double** CrackedVertices_mesh_out, int** CrackedEdgesSize_mesh_out, double** CrackedEdges_mesh_out, int** PreviousNumberingSize_mesh_out, double** PreviousNumbering_mesh_out, int** VerticesSize_geom_out, double** Vertices_geom_out, int** EdgesSize_geom_out, double** Edges_geom_out, int** CornersSize_geom_out, double** Corners_geom_out, int** RequiredVerticesSize_geom_out, double** RequiredVertices_geom_out, int** RequiredEdgesSize_geom_out, double** RequiredEdges_geom_out, int** CrackedEdgesSize_geom_out, double** CrackedEdges_geom_out, int** SubDomainsSize_geom_out, double** SubDomains_geom_out, int* VerticesSize_mesh_in, double* Vertices_mesh_in, int* EdgesSize_mesh_in, double* Edges_mesh_in, int* TrianglesSize_mesh_in, double* Triangles_mesh_in, int* CrackedEdgesSize_mesh_in, double* CrackedEdges_mesh_in, int* VerticesOnGeomEdgeSize_mesh_in, double* VerticesOnGeomEdge_mesh_in, int* VerticesOnGeomVertexSize_mesh_in, double* VerticesOnGeomVertex_mesh_in, int* EdgesOnGeomEdgeSize_mesh_in, double* EdgesOnGeomEdge_mesh_in, int* IssmSegmentsSize_mesh_in, double* IssmSegments_mesh_in, int* VerticesSize_geom_in, double* Vertices_geom_in, int* EdgesSize_geom_in, double* Edges_geom_in, int* CornersSize_geom_in, double* Corners_geom_in, int* RequiredVerticesSize_geom_in, double* RequiredVertices_geom_in, int* RequiredEdgesSize_geom_in, double* RequiredEdges_geom_in, int* CrackedEdgesSize_geom_in, double* CrackedEdges_geom_in, int* SubDomainsSize_geom_in, double* SubDomains_geom_in, double anisomax, double cutoff, double coeff, double errg, double gradation, int Hessiantype, int maxnbv, double maxsubdiv, int Metrictype, int nbjacobi, int nbsmooth, double omega, double power, int verbose, int Crack, int KeepVertices, int splitcorners, double hmin, double hmax, int* hminVerticesSize, double* hminVertices, int* hmaxVerticesSize, double* hmaxVertices, int hVerticesLength, double* hVertices, int* metricSize, double* metric, int* fieldSize, double* field, int* errSize, double* err)

#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  2
#undef NRHS
#define NRHS  3

#endif  /* _BAMG_MESHER_H_ */
