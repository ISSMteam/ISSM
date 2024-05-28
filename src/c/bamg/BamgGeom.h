/*!\file:  BamgGeom.h
 */ 

#ifndef _BAMGGEOM_H_
#define _BAMGGEOM_H_

class BamgGeom{

	public:
		int     VerticesSize[2];
		double* Vertices;
		int     EdgesSize[2];
		double* Edges;
		int     TangentAtEdgesSize[2];
		double* TangentAtEdges;
		int     CornersSize[2];
		double* Corners;
		int     RequiredVerticesSize[2];
		double* RequiredVertices;
		int     RequiredEdgesSize[2];
		double* RequiredEdges;
		int     CrackedEdgesSize[2];
		double* CrackedEdges;
		int     SubDomainsSize[2];
		double* SubDomains;

		BamgGeom();
		~BamgGeom();
};

#endif
