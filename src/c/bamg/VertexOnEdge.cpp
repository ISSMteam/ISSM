#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "VertexOnEdge.h"
#include "Mesh.h"

namespace bamg {

	/*Methods*/
	void VertexOnEdge::Set(const Mesh & Th ,long i,Mesh & ThNew){/*{{{*/
		*this = Th.VertexOnBThEdge[i];  
		v = ThNew.vertices + Th.GetId(v);
	}
	/*}}}*/
	void VertexOnEdge::SetOnBTh(){/*{{{*/
		v->BackgroundEdgeHook=this;
		v->IndexInTriangle=IsVertexOnEdge;  
	}
	/*}}}*/

} 
