#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "GeomVertex.h"
#include "../shared/shared.h"

namespace bamg {

	/*Constructors/Destructors*/
	//See header file

	/*Methods*/
	int  GeomVertex::Corner() const {/*{{{*/
		return type & 4;
	}
	/*}}}*/
	int  GeomVertex::Required()const {/*{{{*/
		// a corner is required
		return type & 6;
	}
	/*}}}*/
	void GeomVertex::SetCorner(){/*{{{*/
		type |= 4;
	}
	/*}}}*/
	void GeomVertex::SetRequired(){/*{{{*/
		type |= 2;
	}
	/*}}}*/

} 
