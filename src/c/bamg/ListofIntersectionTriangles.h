#ifndef _LISTOFINTERSECTIONTRIANGLES_H_
#define _LISTOFINTERSECTIONTRIANGLES_H_

#include "./include.h"

namespace bamg {

	class Triangle;

	class ListofIntersectionTriangles {

		class IntersectionTriangles {

			public: 
				Triangle *t;
				double    bary[3];   // use if t != 0
				R2        x;
				Metric    m;
				double    s;         // curvilinear coordinate
				double    sp;        // length of the previous segment in m
				double    sn;        // length of the next segment in m
		};

		class SegInterpolation {

			public:
				GeomEdge *e;
				double           sBegin  ,sEnd; // abscisse of the seg on edge parameter
				double           lBegin  ,lEnd; // length abscisse set in ListofIntersectionTriangles::Length
				int              last;          // last index in ListofIntersectionTriangles for this Sub seg of edge

				//Methods
				R2 F(double s){ 
					double c01=lEnd-lBegin, c0=(lEnd-s)/c01, c1=(s-lBegin)/c01;
					if (lBegin>s || s>lEnd){
						_error_("lBegin>s || s>lEnd");
					}
					return e->F(sBegin*c0+sEnd*c1);
				}
		};

		public:

			int                    MaxSize;
			int                    Size;
			double                 len;
			int                    state;
			IntersectionTriangles *lIntTria;
			int                    NbSeg;
			int                    MaxNbSeg;
			SegInterpolation      *lSegsI;

			//Constructors/Destructors
			ListofIntersectionTriangles(int n=256,int m=16);
			~ListofIntersectionTriangles();

			//Operators
			IntersectionTriangles & operator[](int i) {return lIntTria[i];}
			operator int&() {return Size;}

			//Methods
			void   Init();
			int    NewItem(Triangle *tt,double d0,double d1,double d2);
			int    NewItem(R2 ,const Metric &);
			void   SplitEdge(Mesh &,const R2 &,const R2 &,int nbegin=0);
			double Length();
			long   NewPoints(BamgVertex *,long &nbv,long maxnbv);
			void   ReShape();
	};

}
#endif
