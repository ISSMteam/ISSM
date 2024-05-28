#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "VertexOnGeom.h"
#include "Mesh.h"
#include "Geometry.h"

namespace bamg {

	/*Constructors/Destructors*/
	VertexOnGeom::VertexOnGeom(){/*{{{*/
		meshvertex=NULL;
		curvilincoord=0;
		gv=0;
	} 
	/*}}}*/
	VertexOnGeom::VertexOnGeom(BamgVertex & m,GeomVertex &g){/*{{{*/
		meshvertex=&m;
		curvilincoord=-1;
		gv=&g;
	}
	/*}}}*/
	VertexOnGeom::VertexOnGeom(BamgVertex & m,GeomEdge &g,double s){/*{{{*/
		meshvertex=&m;
		curvilincoord=s;
		ge=&g;
	}
	/*}}}*/

	/*Methods*/
	void VertexOnGeom::Set(const VertexOnGeom & rec,const Mesh & Th ,Mesh & ThNew){/*{{{*/
		*this = rec;  
		meshvertex = ThNew.vertices + Th.GetId(meshvertex);
		if(gv){
		 if (curvilincoord < 0 )
		  gv = ThNew.Gh.vertices + Th.Gh.GetId(gv);
		 else
		  ge = ThNew.Gh.edges + Th.Gh.GetId(ge);
		}

	}
	/*}}}*/
	int VertexOnGeom::OnGeomVertex()const{/*{{{*/
		return curvilincoord<0;
	}
	/*}}}*/
	int VertexOnGeom::OnGeomEdge() const{/*{{{*/
		return curvilincoord>=0;
	}
	/*}}}*/
	int VertexOnGeom::IsRequiredVertex() {/*{{{*/
		return ((curvilincoord<0 ? (gv?gv->Required():0):0 ));
	}
	/*}}}*/
	void VertexOnGeom::SetOn(){/*{{{*/
		meshvertex->GeomEdgeHook=this;
		meshvertex->IndexInTriangle=IsVertexOnGeom;
	}
	/*}}}*/

} 
