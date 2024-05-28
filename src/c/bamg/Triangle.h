#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include "./include.h"
#include "AdjacentTriangle.h"

namespace bamg {

	class Mesh;
	class BamgVertex;
	class Triangle;

	class Triangle {

		friend class AdjacentTriangle;

		private:
			BamgVertex *vertices[3];        // 3 vertices if t is triangle, t[i] allowed by access function, (*t)[i] if pointer
			Triangle   *adj[3];             // 3 pointers toward the adjacent triangles
			short       AdjEdgeIndex[3];   // edge id in the adjacent triangles. The edge number 1 is the edge number AdjEdgeIndex[1] in the Adjacent triangle 1

		public: 
			long long det; //Integer determinant (twice its area)
			union { 
				Triangle *link;
				long      color;
			};

			//Constructors/Destructors
			Triangle();
			Triangle(Mesh *Th,long i,long j,long k);
			Triangle(BamgVertex *v0,BamgVertex *v1,BamgVertex *v2);

			//Operators
			const BamgVertex & operator[](int i) const {return *vertices[i];};
			BamgVertex & operator[](int i)  {return *vertices[i];};
			const BamgVertex * operator()(int i) const {return vertices[i];};
			BamgVertex * & operator()(int i)  {return vertices[i];};

			//Methods
			void              Echo();
			double            Length() const;
			int               swap(short a1,int=0);
			long              Optim(short a,int =0);
			int               Locked(int a)const;
			int               Hidden(int a)const;
			int               GetAllflag(int a);
			short             NuEdgeTriangleAdj(int i) const;
			AdjacentTriangle  Adj(int i) const;
			Triangle         *TriangleAdj(int i) const;
			void              Renumbering(Triangle   *tb,Triangle *te, long *renu);
			void              Renumbering(BamgVertex *vb,BamgVertex *ve, long *renu);
			void              SetAdjAdj(short a);
			void              SetAdj2(short a,Triangle *t,short aat);
			void              SetSingleVertexToTriangleConnectivity();
			void              SetHidden(int a);
			void              SetLocked(int a);
			void              SetMarkUnSwap(int a);
			void              SetUnMarkUnSwap(int a);

			//Inline methods
			void  Set(const Triangle &,const Mesh &,Mesh &);
			int   In(BamgVertex *v) const { return vertices[0]==v || vertices[1]==v || vertices[2]==v ;}
			BamgVertex* GetVertex(int i){return vertices[i];}; // FIXME: this is used to avoid BamgVertex * operator()

	};

}
#endif
