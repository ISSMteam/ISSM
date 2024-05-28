#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "AdjacentTriangle.h"
#include "Mesh.h"

namespace bamg {

	/*Constructors/Destructors*/
	//See header file

	/*Methods*/
	int  AdjacentTriangle::Locked() const {/*{{{*/
		return t->AdjEdgeIndex[a] & 4;
	}
	/*}}}*/
	int  AdjacentTriangle::GetAllFlag_UnSwap() const {/*{{{*/
		// take all flag except MarkUnSwap
		return t->AdjEdgeIndex[a] & 1012;
	}
	/*}}}*/
	void AdjacentTriangle::SetLock(){/*{{{*/
		t->SetLocked(a);
	}
	/*}}}*/
	AdjacentTriangle AdjacentTriangle::Adj() const {/*{{{*/
		return  t->Adj(a);
	}
	/*}}}*/
	BamgVertex* AdjacentTriangle::EdgeVertex(const int & i) const {/*{{{*/
		return t->vertices[VerticesOfTriangularEdge[a][i]];
	}
	/*}}}*/
	long long & AdjacentTriangle::det() const {/*{{{*/
		return t->det;
	}
	/*}}}*/
	int AdjacentTriangle::swap(){/*{{{*/
		return  t->swap(a);
	}
	/*}}}*/
	void AdjacentTriangle::SetAdj2(const AdjacentTriangle & ta, int l  ){/*{{{*/
		//set Adjacent Triangle of a triangle
		if(t) {
			t->adj[a]=ta.t;
			t->AdjEdgeIndex[a]=ta.a|l;
		}
		if(ta.t) {
			ta.t->adj[ta.a] = t ;
			ta.t->AdjEdgeIndex[ta.a] = a| l ;
		}
	}
	/*}}}*/

} 
