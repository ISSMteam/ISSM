#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "Curve.h"
#include "Geometry.h"

namespace bamg {

	/*Constructors/Destructors*/
	Curve::Curve(){/*{{{*/
		FirstEdge=NULL;
		LastEdge=NULL;
		FirstVertexIndex=0;
		LastVertexIndex=0;
	} 
	/*}}}*/

	/*Methods*/
	void Curve::Set(const Curve & rec,const Geometry & Gh ,Geometry & GhNew){/*{{{*/
		*this = rec;
		FirstEdge = GhNew.edges + Gh.GetId(FirstEdge);    
		LastEdge = GhNew.edges + Gh.GetId(LastEdge); 
	}
	/*}}}*/

} 
