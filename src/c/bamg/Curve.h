#ifndef _CURVE_H_
#define _CURVE_H_

#include "../shared/shared.h"

namespace bamg {

	class GeomEdge;
	class Curve;
	class Geometry;

	class Curve {
		public:
			GeomEdge *FirstEdge; //First edge of the curve
			GeomEdge *LastEdge;  //Last edge of the curve
			int FirstVertexIndex;       //Last vertex index in the last edge
			int LastVertexIndex;        //First Vertex index in the first edge

			//Methods
			Curve();
			void Set(const Curve & rec,const Geometry & Th ,Geometry & ThNew);
	};

}
#endif
