#ifndef _CRACKEDEDGE_H_
#define _CRACKEDEDGE_H_

#include "./typedefs.h"

namespace bamg {

	class Triangle;
	class GeomEdge;
	class Edge;

	class CrackedEdge {

		public:
			Triangle* a;
			Triangle* b; 
			GeomEdge* E;
			Edge* e1;
			Edge* e2;
			double length;
			R2     normal;
			long   indexa[3];
			long   indexb[3];

			//Constructors
			CrackedEdge();
	};

}
#endif
