#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "Edge.h"
#include "Mesh.h"
#include "Geometry.h"
#include "../shared/shared.h"

namespace bamg {

	/*Constructors/Destructors*/

	/*Methods*/
	void Edge::Set(const Mesh & Th ,long i,Mesh & ThNew){ /*{{{*/
		*this = Th.edges[i];
		v[0] = ThNew.vertices + Th.GetId(v[0]);    
		v[1] = ThNew.vertices + Th.GetId(v[1]);
		if (GeomEdgeHook) 
		 GeomEdgeHook =  ThNew.Gh.edges+Th.Gh.GetId(GeomEdgeHook);
		if (adj[0]) adj[0] =   ThNew.edges +   Th.GetId(adj[0]);
		if (adj[1]) adj[1] =   ThNew.edges +   Th.GetId(adj[1]);
	}
	/*}}}*/
	void Edge::Echo(void){ /*{{{*/
		_printf_("Edge:\n");
		_printf_("   pointers towards two vertices: " << v[0] << " " << v[1] << "\n");
		_printf_("   ReferenceNumber = " << ReferenceNumber << "\n");
		_printf_("   GeomEdgeHook = " << GeomEdgeHook << "\n");
		_printf_("   two adjacent edges on the same curve: " << adj[0] << " " << adj[1] << "\n");
	}
	/*}}}*/
	void Edge::Renumbering(BamgVertex *vb,BamgVertex *ve, long *renu){/*{{{*/

		if (v[0] >=vb && v[0] <ve) v[0] = vb + renu[v[0]-vb];
		if (v[1] >=vb && v[1] <ve) v[1] = vb + renu[v[1]-vb];

	}
	/*}}}*/
	int Edge::Intersection(const  Edge & e){ /*{{{*/

		/*some shecks*/
		if (!(adj[0]==&e || adj[1]==&e)){ _error_("Intersection bug"); }
		_assert_(adj[0]==&e || adj[1]==&e);

		return adj[0]==&e?0:1;
	}
	/*}}}*/

} 
