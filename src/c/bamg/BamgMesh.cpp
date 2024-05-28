#include "./bamgobjects.h"
#include "../shared/shared.h"

/*Constructors/Destructors*/
BamgMesh::BamgMesh(){/*{{{*/

	this->VerticesSize[0]=0,                  this->VerticesSize[1]=0;                 this->Vertices=NULL;          this->PreviousNumbering = NULL;
	this->EdgesSize[0]=0,                     this->EdgesSize[1]=0;                    this->Edges=NULL;
	this->TrianglesSize[0]=0,                 this->TrianglesSize[1]=0;                this->Triangles=NULL;

	this->SubDomainsSize[0]=0,                this->SubDomainsSize[1]=0;               this->SubDomains=NULL;
	this->SubDomainsFromGeomSize[0]=0,        this->SubDomainsFromGeomSize[1]=0;       this->SubDomainsFromGeom=NULL;
	this->CrackedVerticesSize[0]=0,           this->CrackedVerticesSize[1]=0;          this->CrackedVertices=NULL;
	this->CrackedEdgesSize[0]=0,              this->CrackedEdgesSize[1]=0;             this->CrackedEdges=NULL;

	this->VerticesOnGeomVertexSize[0]=0,      this->VerticesOnGeomVertexSize[1]=0;     this->VerticesOnGeomVertex=NULL;
	this->VerticesOnGeomEdgeSize[0]=0,        this->VerticesOnGeomEdgeSize[1]=0;       this->VerticesOnGeomEdge=NULL;
	this->EdgesOnGeomEdgeSize[0]=0,           this->EdgesOnGeomEdgeSize[1]=0;          this->EdgesOnGeomEdge=NULL;

	this->IssmEdgesSize[0]=0,                 this->IssmEdgesSize[1]=0;                this->IssmEdges=NULL;
	this->IssmSegmentsSize[0]=0,              this->IssmSegmentsSize[1]=0;             this->IssmSegments=NULL;

	this->ElementConnectivitySize[0]=0,       this->ElementConnectivitySize[1]=0;      this->ElementConnectivity=NULL;
	this->NodalConnectivitySize[0]=0,         this->NodalConnectivitySize[1]=0;        this->NodalConnectivity=NULL;
	this->NodalElementConnectivitySize[0]=0,  this->NodalElementConnectivitySize[1]=0; this->NodalElementConnectivity=NULL;
}
/*}}}*/
BamgMesh::~BamgMesh(){/*{{{*/

	xDelete<double>(this->Vertices);
	xDelete<double>(this->PreviousNumbering);
	xDelete<double>(this->Edges);
	xDelete<double>(this->Triangles);

	xDelete<double>(this->SubDomains);
	xDelete<double>(this->SubDomainsFromGeom);
	xDelete<double>(this->CrackedVertices);
	xDelete<double>(this->CrackedEdges);

	xDelete<double>(this->VerticesOnGeomVertex);
	xDelete<double>(this->VerticesOnGeomEdge);
	xDelete<double>(this->EdgesOnGeomEdge);

	xDelete<double>(this->IssmEdges);
	xDelete<double>(this->IssmSegments);

	xDelete<double>(this->ElementConnectivity);
	xDelete<double>(this->NodalConnectivity);
	xDelete<double>(this->NodalElementConnectivity);
}
/*}}}*/
