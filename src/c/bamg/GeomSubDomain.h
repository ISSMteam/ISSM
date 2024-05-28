#ifndef _GEOMETRICALSUBDOMAIN_H_
#define _GEOMETRICALSUBDOMAIN_H_

#include "./include.h"

namespace bamg {

	class GeomEdge;
	class Geometry;

	class GeomSubDomain {
		public:
			GeomEdge *edge;
			int              direction;   // -1 or 1
			long             ReferenceNumber;

			//Methods
			void Set(const GeomSubDomain &,const Geometry & ,const Geometry &);
	};

}
#endif
