#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "GeomSubDomain.h"
#include "Geometry.h"

namespace bamg {

	/*Constructors/Destructors*/

	/*Methods*/
	void GeomSubDomain::Set(const GeomSubDomain & rec,const Geometry & Gh ,const Geometry & GhNew){/*{{{*/
		*this = rec;
		edge = Gh.GetId(edge) + GhNew.edges;
	}/*}}}*/

} 
