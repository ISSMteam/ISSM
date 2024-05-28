#ifndef _EDGE_H_
#define _EDGE_H_

#include "./BamgVertex.h"
#include "../shared/shared.h"
#include "./GeomEdge.h"

namespace bamg {

	class Mesh;

	class Edge {

		public:
			BamgVertex      *v[2];
			long             ReferenceNumber;
			GeomEdge *GeomEdgeHook;
			Edge            *adj[2];       // the 2 adj edges if on the same curve

			//Operators
			BamgVertex       &operator[](int i){return *v[i];   };
			BamgVertex       *operator()(int     i){return v[i];};
			R2                operator()(double  t) const;// return the point
			const BamgVertex &operator[](int i) const{return *v[i];};

			//Methods
			void Renumbering(BamgVertex *vb,BamgVertex *ve, long *renu);
			int  Intersection(const  Edge & e);
			void Set(const Mesh &,long,Mesh &);
			void Echo(void);

	};
}
#endif
