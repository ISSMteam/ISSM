#ifndef _GEOMETRICALVERTEX_H_
#define _GEOMETRICALVERTEX_H_

#include "./include.h"
#include "BamgVertex.h"

namespace bamg {

	class Geometry;

	class GeomVertex : public BamgVertex { 

		public:
			friend class Geometry;

			int type;

			//Constructors
			GeomVertex():type(0){};

			//Methods
			int  Corner() const;
			int  Required()const;
			void SetCorner();
			void SetRequired();

	};

}
#endif
