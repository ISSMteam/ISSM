#ifndef _GEOMETRICALEDGE_H_
#define _GEOMETRICALEDGE_H_

#include "./include.h"

namespace bamg {

	class GeomVertex;
	class Geometry;

	class GeomEdge {

		public:
			GeomVertex *v[2];
			long               ReferenceNumber;
			long               CurveNumber;
			R2                 tg[2];              // the 2 tangentes (tg[0] =0 => no continuity)
			GeomEdge   *Adj[2];
			int                AdjVertexIndex[2]; // for a given vertex, this gives the index of the vertex in the adjacent edge (0 or 1)
			int                type;

			//Operators
			GeomVertex       & operator[](int i){return *v[i];};
			const GeomVertex & operator[](int i) const { return *v[i];};
			GeomVertex       * operator()(int i){return v[i];};  

			//Methods
			R2     F(double theta) const ; // parametrization of the curve edge
			int    Cracked() const;
			int    TgA()     const;
			int    TgB()     const;
			int    Mark()    const;
			int    Required();
			void   SetCracked();
			void   SetTgA();
			void   SetTgB();
			void   SetMark();
			void   SetUnMark();
			void   SetRequired();
			void   Set(const GeomEdge & rec,const Geometry & Th ,Geometry & ThNew);
	};

}
#endif
