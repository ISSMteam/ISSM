#ifndef _TRIANGLEADJACENT_H_
#define _TRIANGLEADJACENT_H_

#include "./include.h"
#include "./BamgVertex.h"

namespace bamg {

	class Triangle;
	class Triangle;

	class AdjacentTriangle {

		public:
			Triangle* t; //pointer toward the triangle
			int  a;      //Edge number

			//Constructors
			AdjacentTriangle():a(0),t(NULL) {};
			AdjacentTriangle(Triangle* tt,int  aa): t(tt),a(aa &3) {};

			//Operators
			operator Triangle * () const {return t;}
			operator Triangle & () const {return *t;}
			operator int() const {return a;}
			AdjacentTriangle & operator++(){ a= NextEdge[a]; return *this; }
			AdjacentTriangle operator--(){ a= PreviousEdge[a]; return *this; }

			//Methods

			//Methods
			int  Locked() const;
			int  GetAllFlag_UnSwap() const;
			void SetLock();
			void SetAdj2(const AdjacentTriangle &ta, int l=0);
			int  swap();
			AdjacentTriangle Adj() const;
			BamgVertex* EdgeVertex(const int & i) const;
			long long& det() const;
	};
}
#endif
