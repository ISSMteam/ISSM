#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "VertexOnVertex.h"
#include "Mesh.h"

namespace bamg {

	/*Constructors/Destructors*/
	VertexOnVertex::VertexOnVertex() {/*{{{*/
		v=NULL;
		bv=NULL;
	};/*}}}*/
	VertexOnVertex::VertexOnVertex(BamgVertex * w,BamgVertex *bw) :v(w),bv(bw){/*{{{*/

	}/*}}}*/

	/*Methods*/
	void VertexOnVertex::Set(const Mesh &Th ,long i,Mesh &ThNew) { /*{{{*/
		*this = Th.VertexOnBThVertex[i];  
		v     = ThNew.vertices + Th.GetId(v);
	}
	/*}}}*/
	void VertexOnVertex::SetOnBTh(){/*{{{*/
		v->BackgroundVertexHook=bv;v->IndexInTriangle=IsVertexOnVertex;
	}/*}}}*/

} 
