#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "SubDomain.h"
#include "Mesh.h"

namespace bamg {

	/*Constructors/Destructors*/

	/*Methods*/
	void SubDomain::Set(const Mesh & Th ,long i,Mesh & ThNew){/*{{{*/
		*this = Th.subdomains[i];
		if( head-Th.triangles<0 || head-Th.triangles>=Th.nbt){
			_error_("head-Th.triangles<0 || head-Th.triangles>=Th.nbt");
		}
		head = ThNew.triangles + Th.GetId(head) ; 
		if(edge-Th.edges<0 || edge-Th.edges>=Th.nbe){
			_error_("edge-Th.edges<0 || edge-Th.edges>=Th.nbe");
		}
		edge = ThNew.edges+ Th.GetId(edge);
	}
	/*}}}*/

} 
