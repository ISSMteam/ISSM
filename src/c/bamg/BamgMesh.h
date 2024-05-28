/*!\file:  BamgMesh.h
 */ 

#ifndef _BAMGMESH_H_
#define _BAMGMESH_H_

class BamgMesh{

	public:

		int     VerticesSize[2];
		double* Vertices;
		double* PreviousNumbering;
		int     EdgesSize[2];
		double* Edges;
		int     TrianglesSize[2];
		double* Triangles;

		int     VerticesOnGeomVertexSize[2];
		double* VerticesOnGeomVertex;
		int     VerticesOnGeomEdgeSize[2];
		double* VerticesOnGeomEdge;
		int     EdgesOnGeomEdgeSize[2];
		double* EdgesOnGeomEdge;

		int     SubDomainsSize[2];
		double* SubDomains;
		int     SubDomainsFromGeomSize[2];
		double* SubDomainsFromGeom;
		int     CrackedVerticesSize[2];
		double* CrackedVertices;
		int     CrackedEdgesSize[2];
		double* CrackedEdges;

		/*Output for ISSM*/
		int     IssmEdgesSize[2];
		double* IssmEdges;
		int     IssmSegmentsSize[2];
		double* IssmSegments;
		int     ElementConnectivitySize[2];
		double* ElementConnectivity;
		int     NodalConnectivitySize[2];
		double* NodalConnectivity;
		int     NodalElementConnectivitySize[2];
		double* NodalElementConnectivity;

		BamgMesh();
		~BamgMesh();
};

#endif
