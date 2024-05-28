#include "./bamgobjects.h"
#include "../shared/shared.h"

/*Constructors/Destructors*/
BamgGeom::BamgGeom(){/*{{{*/

	this->VerticesSize[0]=0,          this->VerticesSize[1]=0;          this->Vertices=NULL;
	this->EdgesSize[0]=0,             this->EdgesSize[1]=0;             this->Edges=NULL;
	this->TangentAtEdgesSize[0]=0,    this->TangentAtEdgesSize[1]=0;    this->TangentAtEdges=NULL;
	this->CornersSize[0]=0,           this->CornersSize[1]=0;           this->Corners=NULL;
	this->RequiredVerticesSize[0]=0,  this->RequiredVerticesSize[1]=0;  this->RequiredVertices=NULL;
	this->RequiredEdgesSize[0]=0,     this->RequiredEdgesSize[1]=0;     this->RequiredEdges=NULL;
	this->CrackedEdgesSize[0]=0,      this->CrackedEdgesSize[1]=0;      this->CrackedEdges=NULL;
	this->SubDomainsSize[0]=0,        this->SubDomainsSize[1]=0;        this->SubDomains=NULL;

}
/*}}}*/
BamgGeom::~BamgGeom(){/*{{{*/

	xDelete<double>(this->Vertices);
	xDelete<double>(this->Edges);
	xDelete<double>(this->TangentAtEdges);
	xDelete<double>(this->Corners);
	xDelete<double>(this->RequiredVertices);
	xDelete<double>(this->RequiredEdges);
	xDelete<double>(this->CrackedEdges);
	xDelete<double>(this->SubDomains);

}
/*}}}*/
