#ifndef _BAMGVERTEX_H_
#define _BAMGVERTEX_H_

#include "./include.h"
#include "./Metric.h"
#include "./BamgOpts.h"

namespace bamg {

	class Triangle;
	class Mesh;
	class VertexOnGeom;
	class VertexOnEdge;

	class BamgVertex {

		public:

			/*Fields*/
			I2        i;                 // integer coordinates
			R2        r;                 // real coordinates
			Metric    m;
			long      ReferenceNumber;
			long      PreviousNumber;
			short     IndexInTriangle;              // the vertex number in triangle; varies between 0 and 2 in t

			union {
				Triangle     *t;                      // one triangle which is containing the vertex
				long          color;
				BamgVertex   *MeshVertexHook;         // used in geometry BamgVertex to know the Mesh Vertex associated
				VertexOnGeom *GeomEdgeHook;    // if IndexInTriangle == 8; // set with Mesh::SetVertexFieldOn()
				BamgVertex   *BackgroundVertexHook;   // if IndexInTriangle == 16 on Background vertex Mesh::SetVertexFieldOnBTh()
				VertexOnEdge *BackgroundEdgeHook;     // if IndexInTriangle == 32 on Background edge
			};

			/*Operators*/
			operator I2() const {return i;}             // Cast operator
			operator const R2 & () const {return r;}    // Cast operator
			operator Metric () const {return m;}        // Cast operator
			double operator()(R2 x) const { return m.Length(x.x,x.y);} // Get x in the metric m

			/*methods (No constructor and no destructors...)*/
			BamgVertex();
			double Smoothing(Mesh & ,Mesh & ,Triangle  * & ,double =1);
			void   MetricFromHessian(const double Hxx,const double Hyx, const double Hyy, const double smin,const double smax,const double s,const double err,BamgOpts* bamgopts);
			void   Echo();
			int    GetReferenceNumber() const;
			I2     GetIntegerCoordinates() const{return this->i;};// avoid operator I2()
			long   Optim(int =1,int =0); 

			//inline functions
			inline void Set(const BamgVertex &rec,const Mesh & ,Mesh & ){*this=rec;}
	};
}
#endif
